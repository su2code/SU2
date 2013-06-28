/*!
 * \file solution_direct_plasma.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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

#include "../include/solution_structure.hpp"

/*!
 * \class  CPlasmaSolution
 * \brief Initialization of the plasma solution class
 * \author A. Lonkar
 */
CPlasmaSolution::CPlasmaSolution(void) : CSolution() { }

CPlasmaSolution::CPlasmaSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint;
	unsigned short iVar, iDim, iSpecies, iMarker;
	double Vel2 = 0.0;
	string mesh_filename, text_line;
	ifstream restart_file;

	restart = (config->GetRestart() || config->GetRestart_Flow());
	restart_from_Euler = config->GetRestart_Euler2Plasma();
	implicit = config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT;
	centered_scheme = config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED;
	axisymmetric = config->GetAxisymmetric();
	weighted_LS = config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES;


	double Density_Inf_mean, Pressure_Inf_mean, Temperature_Inf_mean, Mach_Inf_mean;
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
	node        = new CVariable*[geometry->GetnPoint()];

	/*--- Define auxiliary vectors for residuals at nodes i & j ---*/
	Residual = new double[nVar];	Residual_Max = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Residual_Chemistry = new double[nVar]; Residual_MomentumExch = new double[nVar];
	Residual_ElecForce = new double[nVar]; Residual_EnergyExch = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar];	Res_Sour = new double[nVar];
	if (axisymmetric) {
		Residual_Axisymmetric = new double[nVar];
	}

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	/*--- Define auxiliary vectors for the solution at nodes i & j ---*/
	Solution = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];

	/*--- Define auxiliary vectors for the geometry ---*/
	Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

	/*--- Define some auxiliary vector related with the undivided lapalacian computation ---*/
	if (centered_scheme) {
		p1_Und_Lapl = new double * [nSpecies];
		p2_Und_Lapl = new double * [nSpecies];
		for (iSpecies =0; iSpecies < nSpecies; iSpecies ++ ) {
			p1_Und_Lapl[iSpecies] = new double [geometry->GetnPoint()];
			p2_Und_Lapl[iSpecies] = new double [geometry->GetnPoint()];
		}
	}

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
		Initialize_Jacobian_Structure(geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs = new double [geometry->GetnPoint()*nVar];
	}

	/*--- Computation of gradients by least squares ---*/
	if (weighted_LS) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar+nSpecies];
		for (iVar = 0; iVar < nVar+nSpecies; iVar++)
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

	CHeatTransfer = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CHeatTransfer[iMarker] = new double [geometry->nVertex[iMarker]];

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
	Gas_Constant_mean    = config->GetGas_Constant();
	SoundSpeed_mean      = sqrt(config->GetGamma()*Gas_Constant_mean*Temperature_Inf_mean);
	Molar_Mass_mean      = config->GetMixtureMolar_Mass();

	/*--- Calculate mean flow density ---*/
	Density_Inf_mean = Pressure_Inf_mean / (Gas_Constant_mean*Temperature_Inf_mean);

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

	/*--- Initialize species quanitites from mean flow properties ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

		/*--- Initialize all species velocites equal to the mean flow velocities ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_Inf[iSpecies][iDim] = Velocity_Inf_mean[iDim];
		}
		Temperature_Inf[iSpecies] = config->GetSpecies_Temperature(iSpecies);
		NumberDensity_ratio       = config->GetInitial_Gas_Composition(iSpecies);
		Molar_Mass                = config->GetMolar_Mass(iSpecies);
		Gas_Constant              = config->GetSpecies_Gas_Constant(iSpecies);
		Gamma                     = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation        = config->GetEnthalpy_Formation(iSpecies);
		NumberDensity             = NumberDensity_ratio*Density_Inf_mean/Molar_Mass_mean;		
		Density_Inf[iSpecies]  	  = Molar_Mass * NumberDensity;		
		Pressure_Inf[iSpecies] = Density_Inf[iSpecies] * Gas_Constant * Temperature_Inf[iSpecies];
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += Velocity_Inf[iSpecies][iDim] * Velocity_Inf[iSpecies][iDim];
		Mach_Inf[iSpecies] 		 = sqrt(Vel2 / (Gamma * Gas_Constant * Temperature_Inf[iSpecies]));
		Energy_vib_Inf[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			CharVibTemp = config->GetCharVibTemp(iSpecies);
			Energy_vib_Inf[iSpecies] = UNIVERSAL_GAS_CONSTANT/Molar_Mass * CharVibTemp / (exp(CharVibTemp/Temperature_Inf[iSpecies]) - 1.0);
		}
		Energy_el_Inf            = 0.0;
		Energy_Inf[iSpecies]     = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*(Gamma-1.0))
																																																			+ 1.0/2.0*Vel2 + Enthalpy_formation + Energy_vib_Inf[iSpecies] + Energy_el_Inf;
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
	int rank = MPI::COMM_WORLD.Get_rank();
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


#ifdef debug	/* Debug start */
	double ionisation = 0.01;
	Temperature_Inf[0] = 810.0;
	Temperature_Inf[1] = 810.0;
	Temperature_Inf[2] = 3900.0;
	double uscalar = 2437.0;
	double M1, M2, M3;
	double R1, R2, R3, Ru;

	M1 = AVOGAD_CONSTANT*config->GetParticle_Mass(0);
	M2 = AVOGAD_CONSTANT*config->GetParticle_Mass(1);
	M3 = AVOGAD_CONSTANT*config->GetParticle_Mass(2);

	Ru = UNIVERSAL_GAS_CONSTANT; R1 = Ru/M1; R2 = Ru/M2; R3 = Ru/M3;
	double *Cv =  new double [nSpecies];

	Gamma = 5.0/3.0;
	Cv[0]  = R1/(Gamma-1.0);
	Cv[1]  = R2/(Gamma-1.0);
	Cv[2]  = R3/(Gamma-1.0);

	Density_Inf[0] = 1.66059E-4;
	Density_Inf[1] = Density_Inf[0]* ionisation *M2/M1;
	Density_Inf[2] = Density_Inf[0]* ionisation *M3/M1;

	Pressure_Inf[0] = Density_Inf[0]*R1*Temperature_Inf[0];
	Pressure_Inf[1] = Density_Inf[1]*R2*Temperature_Inf[1];
	Pressure_Inf[2] = Density_Inf[2]*R3*Temperature_Inf[2];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
		Velocity_Inf[iSpecies][0] = uscalar;
		for (iDim = 1; iDim < nDim; iDim ++)
			Velocity_Inf[iSpecies][iDim] = 0.0;
	}
	double usqr = uscalar*uscalar;

	for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		Energy_Inf[iSpecies] = Cv[iSpecies]*Temperature_Inf[iSpecies] + 0.5*usqr;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
		Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));


	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		cout << " iSpecies = " << iSpecies << endl;
		cout << " Temperature_inf = " << Temperature_Inf[iSpecies] << endl;
		cout << " Pressure_Inf = " << Pressure_Inf[iSpecies] << endl;
		cout << " Density_Inf = " << Density_Inf[iSpecies] << endl;
		cout << " Mach_Inf = " << Mach_Inf[iSpecies] << endl;
		cout << " *********" << endl;
	}

	/* Debug end */

#endif



	/*--- Initialize the solution from configuration file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, nDim, nVar, config);
		}

		/*--- Restart the solution from the specified restart file ---*/
	} else {
		unsigned long index;
		string text_line;

		mesh_filename = config->GetSolution_FlowFileName();
		restart_file.open(mesh_filename.data(), ios::in);

		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}

		/*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
		long *Global2Local;
		Global2Local = new long[geometry->GetGlobal_nPointDomain()];
		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
			Global2Local[iPoint] = -1;
		}
		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
		}

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0;
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
					node[iPoint_Local] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, nDim, nVar, config);
				}
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
			node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, nDim, nVar, config);
		}

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}

}

CPlasmaSolution::~CPlasmaSolution(void) {
	unsigned short iVar, iDim, iMarker, iPoint;

	delete [] Residual;		delete [] Residual_Max;
	delete [] Residual_i;	delete [] Residual_j;
	delete [] Residual_Chemistry;			delete [] Residual_ElecForce;
	delete [] Residual_MomentumExch;	delete[] Residual_EnergyExch;
	delete [] Res_Conv;		delete [] Res_Visc;		delete [] Res_Sour;
	if (axisymmetric)
		delete [] Residual_Axisymmetric;

	delete [] Solution;		delete [] Solution_i; delete [] Solution_j;
	delete [] Vector;			delete [] Vector_i;		delete [] Vector_j;

	if (centered_scheme) {
		for (unsigned short iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
			delete [] p1_Und_Lapl[iSpecies];
			delete [] p2_Und_Lapl[iSpecies];
		}
		delete [] p1_Und_Lapl;
		delete [] p1_Und_Lapl;
	}

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

		delete [] xsol;			delete [] rhs;
	}

	/*--- De-allocate arrays associated with the weighted least-squares gradient calculation method ---*/
	if (weighted_LS) {
		for (iDim = 0; iDim < this->nDim; iDim++)
			delete [] Smatrix[iDim];
		delete [] Smatrix;

		for (iVar = 0; iVar < nVar+nSpecies; iVar++)
			delete [] cvector[iVar];
		delete [] cvector;		
	}

	/*--- De-allocate forces and coefficients on all markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		delete [] CPressure[iMarker];
		delete [] CSkinFriction[iMarker];
		delete [] CHeatTransfer[iMarker];
	}
	delete [] CPressure; delete [] CSkinFriction;	delete [] CHeatTransfer;
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

	/*--- Delete variables stored in at each grid node ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		delete node[iPoint];
	delete [] node;
}

void CPlasmaSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool viscous = (config->GetKind_ViscNumScheme_Plasma());
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {		
		/*--- Compute squared velocity, sound velocity, pressure, and enthalpy ---*/
		node[iPoint]->SetVelocity2();
		node[iPoint]->SetSoundSpeed(config);
		node[iPoint]->SetPressure(config, geometry->node[iPoint]->GetCoord());
		if(node[iPoint]->GetPressure(0) < 0.0 || node[iPoint]->GetPressure(1) < 0.0 || node[iPoint]->GetPressure(2) < 0.0)
			cout << "At iPoint = " <<  iPoint << endl;

		node[iPoint]->SetTemperature_TR(config, geometry->node[iPoint]->GetCoord());
		node[iPoint]->SetTemperature_vib(config, geometry->node[iPoint]->GetCoord());
		node[iPoint]->SetEnthalpy();

		if(viscous) {
			node[iPoint]->SetLaminarViscosity(config);
			node[iPoint]->SetThermalCoeff(config);
		}

		/*--- Initialize the convective and viscous residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResSour_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit)
			node[iPoint]->Set_ResVisc_Zero();
	}
	if (config->GetKind_SourNumScheme_Plasma() !=NONE) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
	}

	/*--- Initialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CPlasmaSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Area, dV, Mean_SoundSpeed, Mean_ProjVel, Lambda, Lambda_1, Lambda_2, Local_Delta_Time, Global_Delta_Time = 1E6;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker, iSpecies, loc;
	double Lambda_iSpecies, Mean_LaminarVisc, Mean_Density;
	double Lambda_Visc_iSpecies,K_v, Local_Delta_Time_Visc;

	K_v = 0.25;
	bool viscous = (config->GetKind_ViscNumScheme_Plasma());
	bool MultipleTimeSteps = (config->MultipleTimeSteps());
	bool centered = (config->GetKind_ConvNumScheme() == SPACE_CENTERED);

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			node[iPoint]->SetMax_Lambda_Inv(0.0, iSpecies);
			node[iPoint]->SetMax_Lambda_Visc(0.0, iSpecies);
			if (centered) node[iPoint]->SetLambda(0.0,iSpecies);
		}

		node[iPoint]->SetVelocity2();
		node[iPoint]->SetPressure(config, geometry->node[iPoint]->GetCoord());
		node[iPoint]->SetSoundSpeed(config);
		if (viscous) {
			node[iPoint]->SetTemperature_TR(config, geometry->node[iPoint]->GetCoord());
			node[iPoint]->SetLaminarViscosity(config);
		}
	}

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

		/*--- Mean Values ---*/
		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal,iSpecies) + node[jPoint]->GetProjVel(Normal,iSpecies));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed(iSpecies) + node[jPoint]->GetSoundSpeed(iSpecies));

			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
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

				Gamma = config->GetSpecies_Gamma(iSpecies);
				Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity(iSpecies) + node[jPoint]->GetLaminarViscosity(iSpecies));
				Mean_Density     = 0.5*(node[iPoint]->GetSolution(loc + 0) + node[jPoint]->GetSolution(loc + 0));

				Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
				Lambda_2 = Gamma*Mean_LaminarVisc/Prandtl_Lam;
				Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

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
			Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {

				/*--- Mean Values ---*/
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal,iSpecies);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed(iSpecies);

				Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
				if (geometry->node[iPoint]->GetDomain()) {

					/*--- Inviscid contribution ---*/
					node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
					if (centered) node[iPoint]->AddLambda(Lambda, iSpecies);

					/*--- Viscous contribution ---*/
					if (viscous) {
						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
						Gamma = config->GetSpecies_Gamma(iSpecies);
						Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity(iSpecies);
						Mean_Density     = node[iPoint]->GetSolution(loc + 0);
						Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
						Lambda_2 = Gamma*Mean_LaminarVisc/Prandtl_Lam;
						Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
						node[iPoint]->AddMax_Lambda_Visc(Lambda, iSpecies);
					}
				}

			}
		}
	}

	if (!MultipleTimeSteps) {
		/*--- Each element uses their own speed ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			dV = geometry->node[iPoint]->GetVolume();
			Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(0);
			for (iSpecies = 1; iSpecies < nSpecies; iSpecies++) {
				Lambda_iSpecies = max(Lambda_iSpecies, node[iPoint]->GetMax_Lambda_Inv(iSpecies));
			}
			Local_Delta_Time = config->GetCFL(iMesh)*dV / Lambda_iSpecies;

			if (viscous) {
				Lambda_Visc_iSpecies = node[iPoint]->GetMax_Lambda_Visc(0);
				for (iSpecies = 1; iSpecies < nSpecies; iSpecies++)
					Lambda_Visc_iSpecies = max(Lambda_Visc_iSpecies, node[iPoint]->GetMax_Lambda_Visc(iSpecies));
				Local_Delta_Time_Visc =  config->GetCFL(iMesh)*K_v*dV*dV / Lambda_Visc_iSpecies;
				Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
			}
			if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
				Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
			}
			else {
				node[iPoint]->SetDelta_Time(Local_Delta_Time);

			}
		}
	}
	else {
		/*--- Each element uses their own speed, steady state simulation ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			dV = geometry->node[iPoint]->GetVolume();

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(iSpecies);
				Local_Delta_Time = config->GetCFL(iMesh)*dV / Lambda_iSpecies;
				if (viscous) {
					Lambda_Visc_iSpecies = node[iPoint]->GetMax_Lambda_Visc(iSpecies);
					Local_Delta_Time_Visc =  config->GetCFL(iMesh)*K_v*dV*dV / Lambda_Visc_iSpecies;
					Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
				}
				node[iPoint]->SetDelta_Time(Local_Delta_Time,iSpecies);
			}
		}
	}

	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}
}

void CPlasmaSolution::Centered_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) { 	unsigned long iEdge, iPoint, jPoint;

		bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
		bool dissipation = ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit);
		bool high_order_diss = ((config->GetKind_Centered() == JST) && (iMesh == MESH_0));
		unsigned short iSpecies;

		/*--- Artificial dissipation preprocessing ---*/
		if (dissipation && high_order_diss) {
			SetDissipation_Switch(geometry, solution_container, config);
			SetUndivided_Laplacian(geometry, config);
		}

		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			/*--- Points in edge, set normal vectors, and number of neighbors ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());
			solver->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());

			/*--- Set conservative variables w/o reconstruction ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				solver->SetPressure(node[iPoint]->GetPressure(iSpecies), node[jPoint]->GetPressure(iSpecies), iSpecies);
				solver->SetSoundSpeed(node[iPoint]->GetSoundSpeed(iSpecies), node[jPoint]->GetSoundSpeed(iSpecies),iSpecies);
				solver->SetEnthalpy(node[iPoint]->GetEnthalpy(iSpecies), node[jPoint]->GetEnthalpy(iSpecies),iSpecies);
				solver->SetLambda(node[iPoint]->GetLambda(iSpecies), node[jPoint]->GetLambda(iSpecies),iSpecies);
			}


			if (dissipation && high_order_diss) {
				solver->SetUndivided_Laplacian(node[iPoint]->GetUnd_Lapl(), node[jPoint]->GetUnd_Lapl());
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					solver->SetSensor(node[iPoint]->GetSensor(iSpecies), node[jPoint]->GetSensor(iSpecies), iSpecies);
				}
			}


			/*--- Compute residuals ---*/
			solver->SetResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, config);

			/*--- Update convective and artificial dissipation residuals ---*/
			node[iPoint]->AddRes_Conv(Res_Conv);
			node[jPoint]->SubtractRes_Conv(Res_Conv);
			if (dissipation) {
				node[iPoint]->AddRes_Visc(Res_Visc);
				node[jPoint]->SubtractRes_Visc(Res_Visc);
			}

			/*--- Set implicit stuff ---*/
			if (implicit) {
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
				Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
				Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);

			}
		}
}


void CPlasmaSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool high_order_diss = (( (config->GetKind_Upwind() == ROE_2ND) || (config->GetKind_Upwind() == SW_2ND))&& (iMesh == MESH_0));
	double *mixed_U_i, *mixed_U_j;
	//SetSolution_Gradient_GG(geometry);

	mixed_U_i = new double[nVar];
	mixed_U_j = new double[nVar];

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
	}

	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Set conservative variables w/o reconstruction ---*/
		U_i = node[iPoint]->GetSolution();
		U_j = node[jPoint]->GetSolution();

		solver->SetConservative(U_i, U_j);

		if (high_order_diss) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}

			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_j = node[jPoint]->GetLimiter(); Limiter_i = node[iPoint]->GetLimiter();

        /*--- If the node belongs to the boundary, revert to first order ---*/
        if (geometry->node[iPoint]->GetBoundary() || geometry->node[jPoint]->GetBoundary())
          for (iVar = 0; iVar < nVar; iVar++) {
            Limiter_i[iVar] = 0.0;
            Limiter_j[iVar] = 0.0;
          }
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
//          Solution_i[iVar] = U_i[iVar] + Project_Grad_i*Limiter_i[iVar]*1E-1;
//					Solution_j[iVar] = U_j[iVar] + Project_Grad_j*Limiter_j[iVar]*1E-1;
				}
			}

			/*-- New implementation ---*/
			/*      double alpha, Pressure_i, Pressure_j, gradP, w;
      unsigned short iSpecies, loc;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
        else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
        alpha = 0.5;
        Pressure_i = node[iPoint]->GetPressure(iSpecies);
        Pressure_j = node[jPoint]->GetPressure(iSpecies);          
        gradP = fabs(Pressure_j - Pressure_i) / min(Pressure_j, Pressure_i);
        w = 1.0 / (alpha*gradP*alpha*gradP + 1.0);
        if (iEdge == 16) {
          cout << "w: " << w << endl;
          cin.get();
        }
        for (iVar = 0; iVar < nDim+2; iVar++) {
          Solution_i[loc+iVar] -= (1-w)*Project_Grad_i*Limiter_i[loc+iVar];          
          Solution_j[loc+iVar] -= (1-w)*Project_Grad_j*Limiter_j[loc+iVar];          
        }
        if (iSpecies < nDiatomics) {
          Solution_i[loc+nDim+2] -= (1-w)*Project_Grad_i*Limiter_i[loc+nDim+2];          
          Solution_j[loc+nDim+2] -= (1-w)*Project_Grad_j*Limiter_j[loc+nDim+2];           
        }          
      }*/

			/*--- Set conservative variables with reconstruction ---*/
			for (iVar = 0; iVar < nVar; iVar ++) {
				mixed_U_i[iVar] = 0.0*U_i[iVar] + 1.0*Solution_i[iVar];
				mixed_U_j[iVar] = 0.0*U_j[iVar] + 1.0*Solution_j[iVar];
			}
			solver->SetConservative(mixed_U_i, mixed_U_j);
		}

		/*--- Compute the residual ---*/
		solver->SetResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

		/*--- Update residual value ---*/
		node[iPoint]->AddRes_Conv(Res_Conv);
		node[jPoint]->SubtractRes_Conv(Res_Conv);

		/*--- Set implicit jacobians ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);

		}
	}
	delete [] mixed_U_i; delete [] mixed_U_j;
}

void CPlasmaSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge, iSpecies;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		/*--- Compute gradients (not in the preprocessing, because we want to be sure that we have
		 the gradient of the primitive varibles, not the conservative) ---*/

		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimVar_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimVar_Gradient_LS(geometry, config);

		for (iEdge = 0; iEdge <  geometry->GetnEdge(); iEdge++) {

			/*--- Points, coordinates and normal vector in edge ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());

			/*--- Conservative variables, and primitive Variables w/o reconstruction ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),node[jPoint]->GetGradient_Primitive());

			/*--- Laminar and eddy viscosity ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				solver->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(iSpecies), node[jPoint]->GetLaminarViscosity(iSpecies), iSpecies);
				solver->SetEddyViscosity(node[iPoint]->GetEddyViscosity(iSpecies), node[jPoint]->GetEddyViscosity(iSpecies), iSpecies);
			}
			/*--- Compute and update residual ---*/
			solver->SetResidual(Res_Visc, Jacobian_i, Jacobian_j, config);

			node[iPoint]->SubtractRes_Visc(Res_Visc);
			node[jPoint]->AddRes_Visc(Res_Visc);

			/*--- Implicit part ---*/
			if (implicit) {
				Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
				Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
				Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
				Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
			}
		}
	}
}

void CPlasmaSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {
	unsigned short iVar, jVar, iSpecies;
	unsigned long iPoint;
	double *Efield;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool axisymmetric = config->GetAxisymmetric();

	double *Temperature_tr_i, *Temperature_vib_i, *Pressure_i;

	Temperature_tr_i = new double [nSpecies];
	Temperature_vib_i = new double [nSpecies];
	Pressure_i = new double[nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Temperature_tr_i[iSpecies] = 0.0;
		Temperature_vib_i[iSpecies] = 0.0;
		Pressure_i[iSpecies] = 0.0;
	}

	Efield = new double [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim ++)
		Efield[iDim] = 0.0;

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

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

		if (config->GetElectricSolver())
			Efield = node[iPoint]->GetElectricField();

		solver->SetElectricField(Efield);

		/*--- Set y coordinate ---*/
		solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());

		/*--- Set solution  ---*/
		solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

		/*--- Set control volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Set temperature ---*/
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
			Temperature_tr_i[iSpecies] = node[iPoint]->GetTemperature_TR(iSpecies);
			Temperature_vib_i[iSpecies] = node[iPoint]->GetTemperature_vib(iSpecies);
			Pressure_i[iSpecies] = node[iPoint]->GetPressure(iSpecies);
		}
		solver->SetTemperature_TR(Temperature_tr_i, Temperature_tr_i);
		solver->SetTemperature_vib(Temperature_vib_i, Temperature_vib_i);

		/*--- Set pressure ---*/
		solver->SetPressure(Pressure_i, Pressure_i);

		/*--- Compute Residual ---*/
		if (config->GetKind_GasModel() == ARGON || config->GetKind_GasModel() == AIR21) {
			solver->SetResidual(Residual, Jacobian_i, config);
			node[iPoint]->SubtractRes_Sour(Residual);
			if (implicit)
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
		}

		else if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2 
				|| config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {

			/*--- Axisymmetric source terms ---*/
			if (axisymmetric) {
				solver->SetResidual_Axisymmetric(Residual_Axisymmetric, config);
				solver->SetJacobian_Axisymmetric(Jacobian_Axisymmetric, config);
				node[iPoint]->AddRes_Sour(Residual_Axisymmetric);
				if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_Axisymmetric);	
			}

			/*--- Chemistry source terms ---*/
			solver->SetResidual_Chemistry(Residual_Chemistry, config);
			solver->SetJacobian_Chemistry(Jacobian_Chemistry, config);
			node[iPoint]->SubtractRes_Sour(Residual_Chemistry);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_Chemistry);			

/*      if (iPoint == 17) {
        cout << "U_i: " << endl;
        for (iVar = 0; iVar < nVar; iVar++)
          cout << Solution_i[iVar] << endl;
        cout << endl << endl;
        cout << "Residual Chemistry: " << endl;
        for (iVar = 0; iVar < nVar; iVar++)
          cout << Residual_Chemistry[iVar] << endl;
        cout << endl << endl << "Jacobian Chemistry: " << endl;
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            cout << Jacobian_Chemistry[iVar][jVar] << "\t";
          }
          cout << endl;
        }
        cin.get();
      }  */

			/*--- Electric force source terms ---*/
			/*		solver->SetResidual_ElecForce(Residual_ElecForce, config);
			 solver->SetJacobian_ElecForce(Jacobian_ElecForce, config);
			 node[iPoint]->SubtractRes_Sour(Residual_ElecForce);
			 if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ElecForce);*/

			/*--- Momentum exchange source terms ---*/
			solver->SetResidual_MomentumExch(Residual_MomentumExch, config);
			solver->SetJacobian_MomentumExch(Jacobian_MomentumExch, config);			
			node[iPoint]->SubtractRes_Sour(Residual_MomentumExch);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_MomentumExch);

			/*--- Energy exchange source terms ---*/
			solver->SetResidual_EnergyExch(Residual_EnergyExch, config);
			solver->SetJacobian_EnergyExch(Jacobian_EnergyExch, config);
			node[iPoint]->SubtractRes_Sour(Residual_EnergyExch);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_EnergyExch);
		}
	}
	delete[] Temperature_tr_i;
	delete[] Temperature_vib_i;
	delete[] Pressure_i;
	delete [] Efield;

}

/*!
 * \method Copy_Zone_Solution
 * \brief Copy solution from solver 1 into solver 2
 * \author A. Lonkar
 */
void CPlasmaSolution::Copy_Zone_Solution(CSolution ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config,
		CSolution ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {
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


void CPlasmaSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

	unsigned long iPoint, jPoint, iEdge, Point_Normal = 0, iVertex;
	double Pressure_i = 0, Pressure_j = 0, *Normal;
	unsigned short iVar, iMarker, iSpecies, loc, nVar_Species;
	double *Diff = new double[nVar];
	double *U_halo = new double[nVar];

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
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

					Point_Normal = geometry->vertex[iMarker][iVertex]->GetClosest_Neighbor();

					/*--- Interpolate & compute difference in the conserved variables ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL
								|| config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL) {
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

void CPlasmaSolution::Inviscid_Forces(CGeometry *geometry, CConfig *config) {

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
	Total_CT = 0.0; Total_CQ = 0.0;
	AllBound_CDrag_Inv = 0.0; AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0;
	AllBound_CMx_Inv = 0.0; AllBound_CMy_Inv = 0.0; AllBound_CMz_Inv = 0.0;
	AllBound_CFx_Inv = 0.0; AllBound_CFy_Inv = 0.0; AllBound_CFz_Inv = 0.0;
	AllBound_CEff_Inv = 0.0; AllBound_CMerit_Inv = 0.0;
	AllBound_CNearFieldOF_Inv = 0.0;
	AllBound_CT_Inv = 0.0; AllBound_CQ_Inv = 0.0;

	/*--- Loop over the Euler and Navier-Stokes markers ---*/
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary   = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);

		if ((Boundary == EULER_WALL) || (Boundary == NO_SLIP_WALL) || (Boundary == NEARFIELD_BOUNDARY)) {

			for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;

			MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;

			NFPressOF = 0.0; PressInviscid = 0.0;

			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

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

void CPlasmaSolution::Viscous_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short Boundary, iMarker, iDim, jDim, iSpecies, loc;
	double **Tau, Delta = 0.0, Viscosity, **Grad_PrimVar, div_vel, *Normal, **TauElem;
	double  *Kappa, *TauTangent, Area, WallShearStress, *TauNormal;
	//double factor, Wall_heat_transfer;
	//double cp, Laminar_Viscosity,  heat_flux_factor, GradTemperature;
	double RefDensity, RefVel2;
	double CSkinFriction_iSpecies;// CHeatTransfer_iSpecies;

	//	double Gas_Constant = config->GetGas_Constant();
	RefDensity = Density_Inf[0];


	//cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

	/*--- Vector and variables initialization ---*/
	Kappa      = new double [nDim];
	TauElem    = new double* [nSpecies];
	TauTangent = new double [nDim];
	Tau        = new double* [nDim];
	TauNormal = new double [nSpecies];

	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim]   = new double [nDim];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
		TauElem[iSpecies] = new double [nDim];


	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0 ; jDim < nDim; jDim++) {
			Delta = 0.0; if (iDim == jDim) Delta = 1.0;
		}

	/*--- Loop over the Navier-Stokes markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);

		if (Boundary == NO_SLIP_WALL) {

			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {

					Grad_PrimVar = node[iPoint]->GetGradient_Primitive();

					/*--- Compute viscous foces ---*/
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
					for (iDim = 0; iDim < nDim; iDim++) Kappa[iDim] = Normal[iDim]/Area;

					CSkinFriction[iMarker][iVertex] = 0.0;
					//CHeatTransfer[iMarker][iVertex] = 0.0;

					for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
						Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

						div_vel = 0.0;
						for (iDim = 0; iDim < nDim; iDim++)
							div_vel += Grad_PrimVar[loc + iDim+1][iDim];

						for (iDim = 0; iDim < nDim; iDim++)
							for (jDim = 0 ; jDim < nDim; jDim++) {
								Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[loc + jDim+1][iDim] +
										Grad_PrimVar[loc + iDim+1][jDim]) - TWO3*Viscosity*div_vel*Delta;
							}

						for (iDim = 0; iDim < nDim; iDim++) {
							TauElem[iSpecies][iDim] = 0.0;
							for (jDim = 0; jDim < nDim; jDim++)
								TauElem[iSpecies][iDim] += Tau[iDim][jDim]*Kappa[jDim] ;
						}

						/*--- Compute wall shear stress, and skin friction coefficient ---*/
						TauNormal[iSpecies] = 0.0;
						for (iDim = 0; iDim < nDim; iDim++)
							TauNormal[iSpecies] += TauElem[iSpecies][iDim] * Kappa[iDim];

						for (iDim = 0; iDim < nDim; iDim++)
							TauTangent[iDim] = TauElem[iSpecies][iDim] - TauNormal[iSpecies] * Kappa[iDim];

						WallShearStress = 0.0;
						for (iDim = 0; iDim < nDim; iDim++)
							WallShearStress += TauTangent[iDim]*TauTangent[iDim];

						RefDensity		= Density_Inf[iSpecies];
						RefVel2 = 0.0;
						for (iDim = 0; iDim < nDim; iDim ++)
							RefVel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];

						CSkinFriction_iSpecies = sqrt(WallShearStress) / (0.5*RefDensity*RefVel2);
						CSkinFriction[iMarker][iVertex] += CSkinFriction_iSpecies;

						/*	Laminar_Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies) ;
						gas_Constant = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));
						cp = (Gamma / Gamma_Minus_One) * gas_Constant;
						heat_flux_factor = cp * Laminar_Viscosity/PRANDTL;
						GradTemperature = 0.0;
						for (iDim = 0; iDim < nDim; iDim++)
							GradTemperature +=  Grad_PrimVar[loc + 0][iDim]*(-Normal[iDim]);

						Wall_heat_transfer = heat_flux_factor*GradTemperature;
						CHeatTransfer_iSpecies = Wall_heat_transfer/(0.5*RefDensity*RefVel2);
						CHeatTransfer[iMarker][iVertex] += CHeatTransfer_iSpecies;
						 */
					}
				}
			}
		}
	}

	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
	delete [] Kappa;
	delete [] TauTangent;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
		delete [] TauElem[iSpecies];

	delete [] TauElem;

}


void CPlasmaSolution::SetDissipation_Switch(CGeometry *geometry, CSolution **solution_container, CConfig *config) {

	unsigned long iEdge, iPoint, jPoint;
	double Pressure_i, Pressure_j;
	unsigned short iSpecies;
	/*--- Reset variables to store the undivided pressure ---*/
	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++ ) {
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
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
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
			node[iPoint]->SetSensor(fabs(p1_Und_Lapl[iSpecies][iPoint]) / p2_Und_Lapl[iSpecies][iPoint], iSpecies);

		}
	}
}

void CPlasmaSolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container,
		CConfig *config, unsigned short iRKStep) {
	double *Residual, *Res_TruncError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;
		Res_TruncError = node[iPoint]->GetRes_TruncError();
		Residual = node[iPoint]->GetResidual();
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -(Residual[iVar]+Res_TruncError[iVar])*Delta*RK_AlphaCoeff);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
		}
	}

#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
#endif

}

void CPlasmaSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *Residual, *Res_TruncError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;

	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;
		Res_TruncError = node[iPoint]->GetRes_TruncError();
		Residual = node[iPoint]->GetResidual();

		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -(Residual[iVar])*Delta);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
		}
	}

#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
#endif

}

void CPlasmaSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar, iSpecies;
	unsigned long iPoint, total_index;
	double Delta, Res, *local_ResConv, *local_ResVisc, *local_ResSour, *local_Res_TruncError, Vol;
	double *Species_Delta;
	bool MultipleTimeSteps = (config->MultipleTimeSteps());

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	Species_Delta = new double [nSpecies];

	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		local_Res_TruncError = node[iPoint]->GetRes_TruncError();
		local_ResConv = node[iPoint]->GetResConv();
		local_ResVisc = node[iPoint]->GetResVisc();
		local_ResSour = node[iPoint]->GetResSour();
		Vol = geometry->node[iPoint]->GetVolume();

		if (!MultipleTimeSteps) {
			/*--- Modify matrix diagonal to assure diagonal dominance ---*/
			Delta = geometry->node[iPoint]->GetVolume()/node[iPoint]->GetDelta_Time();
			Jacobian.AddVal2Diag(iPoint,Delta);
		} else {
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				Species_Delta[iSpecies] = geometry->node[iPoint]->GetVolume()/node[iPoint]->GetDelta_Time(iSpecies);
			}
			Jacobian.AddVal2Diag(iPoint, Species_Delta, nDim);
		}

		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
			Res = local_ResConv[iVar]+local_ResVisc[iVar]+local_ResSour[iVar];
			rhs[total_index] = -(Res+local_Res_TruncError[iVar]);
			xsol[total_index] = 0.0;
			AddRes_Max( iVar, Res*Res*Vol );
		}

	}

	/*--- Solve the system ---*/
	if(config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2 || 
			config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {

		/*--- Solve the linear system (Stationary iterative methods) ---*/
		if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) 
			Jacobian.SGSSolution(rhs, xsol, config->GetLinear_Solver_Error(), 
					config->GetLinear_Solver_Iter(), false, geometry, config);

		if (config->GetKind_Linear_Solver() == LU_SGS) 
			Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);

		/*--- Solve the linear system (Krylov subspace methods) ---*/
		if ((config->GetKind_Linear_Solver() == BCGSTAB) || 
				(config->GetKind_Linear_Solver() == GMRES)) {

			CSysVector rhs_vec(geometry->GetnPoint(), geometry->GetnPointDomain(), nVar, rhs);
			CSysVector sol_vec(geometry->GetnPoint(), geometry->GetnPointDomain(), nVar, xsol);

			CMatrixVectorProduct* mat_vec = new CSparseMatrixVectorProduct(Jacobian);
			CSolutionSendReceive* sol_mpi = new CSparseMatrixSolMPI(Jacobian, geometry, config);

			CPreconditioner* precond = NULL;
			if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
				Jacobian.BuildJacobiPreconditioner();
				precond = new CJacobiPreconditioner(Jacobian);			
			}
			else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
				Jacobian.BuildJacobiPreconditioner();
				precond = new CLineletPreconditioner(Jacobian);
			}
			else if (config->GetKind_Linear_Solver_Prec() == NO_PREC) 
				precond = new CIdentityPreconditioner();

			CSysSolve system;
			if (config->GetKind_Linear_Solver() == BCGSTAB)
				system.BCGSTAB(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(), 
						config->GetLinear_Solver_Iter(), false);
			else if (config->GetKind_Linear_Solver() == GMRES)
				system.FlexibleGMRES(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(), 
						config->GetLinear_Solver_Iter(), false);

			sol_vec.CopyToArray(xsol);
			delete mat_vec; 
			delete precond;
			delete sol_mpi;
		}

		/*--- Update solution (system written in terms of increments) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);

	}

	if(config->GetKind_GasModel() == AIR21) {
		if (((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) ||
				(config->GetUnsteady_Simulation() == TIME_STEPPING)) {
			/*--- If it is a unsteady problem, the linear system must be converged ---*/

			Jacobian.SGSSolution(rhs, xsol, 1e-9, 99999, true, geometry, config);
		}
		else 	Jacobian.SGSSolution(rhs, xsol, 1e-9, 500, true, geometry, config);


		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);
	}

	if(config->GetKind_GasModel() == ARGON) {
		if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) Jacobian.SGSSolution(rhs, xsol, 1e-9, 500, true, geometry, config);
		if (config->GetKind_Linear_Solver() == LU_SGS) Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);

		/*--- Solve the linear system (Krylov subspace methods) ---*/
		if ((config->GetKind_Linear_Solver() == BCGSTAB) || 
				(config->GetKind_Linear_Solver() == GMRES)) {

			CSysVector rhs_vec(geometry->GetnPoint(), geometry->GetnPointDomain(), nVar, rhs);
			CSysVector sol_vec(geometry->GetnPoint(), geometry->GetnPointDomain(), nVar, xsol);

			CMatrixVectorProduct* mat_vec = new CSparseMatrixVectorProduct(Jacobian);
			CSolutionSendReceive* sol_mpi = new CSparseMatrixSolMPI(Jacobian, geometry, config);

			CPreconditioner* precond = NULL;
			if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
				Jacobian.BuildJacobiPreconditioner();
				precond = new CJacobiPreconditioner(Jacobian);			
			}
			else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
				Jacobian.BuildJacobiPreconditioner();
				precond = new CLineletPreconditioner(Jacobian);
			}
			else if (config->GetKind_Linear_Solver_Prec() == NO_PREC) 
				precond = new CIdentityPreconditioner();

			CSysSolve system;
			if (config->GetKind_Linear_Solver() == BCGSTAB)
				system.BCGSTAB(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, 1e-9, 500, false);
			else if (config->GetKind_Linear_Solver() == GMRES)
				system.FlexibleGMRES(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, 1e-9, 500, false);		

			sol_vec.CopyToArray(xsol);
			delete mat_vec; 
			delete precond;
			delete sol_mpi;
		}


		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);
		}

#ifdef FirstProblem
		unsigned short interior = 0;

		if (nDim ==2 ) {
			interior = 10;
			/*--- Update solution (system written in terms of increments) ---*/

			for (iPoint = interior; iPoint < geometry->GetnPointDomain(); iPoint++)
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);

			for (iPoint = 0; iPoint < interior; iPoint++)
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->AddSolution(iVar,xsol[(interior+iPoint)*nVar+iVar]);


			for (iPoint = 0; iPoint < interior; iPoint++)
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->SetSolution(iVar,node[iPoint+interior]->GetSolution(iVar));

		}

		if (nDim ==3 ) {
			unsigned short Start = 400;
			unsigned short iJump  = 1;
			unsigned long iStart;
			/*--- Update solution (system written in terms of increments) ---*/

			for (iPoint = interior; iPoint < geometry->GetnPointDomain(); iPoint++)
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);


			for (iPoint = 0; iPoint < interior; iPoint++)
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->AddSolution(iVar,xsol[(interior+iPoint)*nVar+iVar]);


			for (iJump = 0; iJump <5; iJump ++) {
				iStart = Start + iJump * 405;
				for (iPoint = iStart; iPoint < iStart + 5; iPoint++)
					for (iVar = 0; iVar < nVar; iVar++)
						node[iPoint]->SetSolution(iVar,node[iPoint-5]->GetSolution(iVar));
			}
		}

# endif
	}

#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
#endif

}

void CPlasmaSolution::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iVar, iMarker, iSpecies, loc;
	double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
	Partial_Gradient, Partial_Res, *Normal;

	PrimVar_Vertex = new double [nVar+nSpecies];
	PrimVar_i = new double [nVar+nSpecies];
	PrimVar_j = new double [nVar+nSpecies];

	/*--- Set Gradient_Primitive to Zero ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero();

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+4)*iSpecies;
			else loc = (nDim+4)*nDiatomics + (nDim+3)*(iSpecies-nDiatomics);

			PrimVar_i[loc + 0] = node[iPoint]->GetTemperature(iSpecies);
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_i[loc + iDim+1] = node[iPoint]->GetVelocity(iDim,iSpecies);
			PrimVar_i[loc + nDim+1] = node[iPoint]->GetPressure(iSpecies);
			PrimVar_i[loc + nDim+2] = node[iPoint]->GetDensity(iSpecies);

			PrimVar_j[loc + 0] = node[jPoint]->GetTemperature(iSpecies);
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_j[loc + iDim+1] = node[jPoint]->GetVelocity(iDim,iSpecies);
			PrimVar_j[loc + nDim+1] = node[jPoint]->GetPressure(iSpecies);
			PrimVar_j[loc + nDim+2] = node[jPoint]->GetDensity(iSpecies);
		}

		Normal = geometry->edge[iEdge]->GetNormal();
		for (iVar = 0; iVar < nVar+nSpecies; iVar++) {
			PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
			for (iDim = 0; iDim < nDim; iDim++) {
				Partial_Res = PrimVar_Average*Normal[iDim];
				if (geometry->node[iPoint]->GetDomain())
					node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
				if (geometry->node[jPoint]->GetDomain())
					node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
			}
		}
	}

	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+4)*iSpecies;
					else loc = (nDim+4)*nDiatomics + (nDim+3)*(iSpecies-nDiatomics);

					PrimVar_Vertex[loc + 0] = node[iPoint]->GetTemperature(iSpecies);
					for (iDim = 0; iDim < nDim; iDim++)
						PrimVar_Vertex[loc + iDim+1] = node[iPoint]->GetVelocity(iDim,iSpecies);
					PrimVar_Vertex[loc + nDim+1] = node[iPoint]->GetPressure(iSpecies);
					PrimVar_Vertex[loc + nDim+2] = node[iPoint]->GetDensity(iSpecies);
				}

				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				for (iVar = 0; iVar < nVar+nSpecies; iVar++)
					for (iDim = 0; iDim < nDim; iDim++) {
						Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
						node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
					}
			}
		}
	}

	/*--- Update gradient value ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar+nSpecies; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar,iDim) / (geometry->node[iPoint]->GetVolume()+EPS);
				//				if (Partial_Gradient > config->GetPrimGrad_Threshold()) Partial_Gradient = config->GetPrimGrad_Threshold();
				//				if (Partial_Gradient < -config->GetPrimGrad_Threshold()) Partial_Gradient = -config->GetPrimGrad_Threshold();
				node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
			}
		}
	}


	delete [] PrimVar_Vertex;
	delete [] PrimVar_i;
	delete [] PrimVar_j;
}

void CPlasmaSolution::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, jDim, iNeigh, iSpecies, loc = 0;
	unsigned long iPoint, jPoint;
	double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, weight, product;

	PrimVar_i = new double [nVar+nSpecies];
	PrimVar_j = new double [nVar+nSpecies];

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Coord_i = geometry->node[iPoint]->GetCoord();

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+4)*iSpecies;
			else loc = (nDim+4)*nDiatomics + (nDim+3)*(iSpecies-nDiatomics);

			PrimVar_i[loc + 0] = node[iPoint]->GetTemperature(iSpecies);
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_i[loc + iDim+1] = node[iPoint]->GetVelocity(iDim,iSpecies);
			PrimVar_i[loc + nDim+1] = node[iPoint]->GetPressure(iSpecies);
			PrimVar_i[loc + nDim+2] = node[iPoint]->GetDensity(iSpecies);
		}

		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nVar+nSpecies; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iVar][iDim] = 0.0;
		r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+4)*iSpecies;
				else loc = (nDim+4)*nDiatomics + (nDim+3)*(iSpecies-nDiatomics);

				PrimVar_j[loc + 0] = node[jPoint]->GetTemperature(iSpecies);
				for (iDim = 0; iDim < nDim; iDim++)
					PrimVar_j[loc + iDim+1] = node[jPoint]->GetVelocity(iDim,iSpecies);
				PrimVar_j[loc + nDim+1] = node[jPoint]->GetPressure(iSpecies);
				PrimVar_j[loc + nDim+2] = node[jPoint]->GetDensity(iSpecies);
			}
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
			for (iVar = 0; iVar < nVar+nSpecies; iVar++)
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
		for (iVar = 0; iVar < nVar+nSpecies; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				product = 0.0;
				for (jDim = 0; jDim < nDim; jDim++) {
					product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
				}
				node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
			}

		}

	}

	delete [] PrimVar_i;
	delete [] PrimVar_j;
}

void CPlasmaSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iSpecies, iVar, loc;
	double Pressure, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {

				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

				Pressure = node[iPoint]->GetPressure(iSpecies);

				Residual[loc+0] = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[loc + iDim+1] =  Pressure*UnitaryNormal[iDim]*Area;
				Residual[loc+nDim+1] = 0.0;
				if ( iSpecies < nDiatomics )
					Residual[loc+nDim+2] = 0.0;
			}

			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);

			/*--- In case we are doing a implicit computation ---*/

			if (implicit) {
				double a2;
				double phi;
        double Energy_el;
        Energy_el = 0.0;
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
					a2 = config->GetSpecies_Gamma(iSpecies)-1.0;
					phi = a2*(0.5*node[iPoint]->GetVelocity2(iSpecies) - config->GetEnthalpy_Formation(iSpecies) - Energy_el);

					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[loc + 0][iVar] = 0.0;
						Jacobian_i[loc + nDim+1][iVar] = 0.0;
					}
          if (iSpecies < nDiatomics)
            for (iVar = 0; iVar < nVar; iVar++)
              Jacobian_i[loc+nDim+2][iVar] = 0.0;

					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
            if (iSpecies < nDiatomics)
              Jacobian_i[loc+iDim+1][loc+nDim+2] = a2*Normal[iDim];
					}
				}
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

void CPlasmaSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
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
					node[iPoint]->SetVel_ResConv_Zero(iSpecies);
					node[iPoint]->SetVel_ResVisc_Zero(iSpecies);
					node[iPoint]->SetVelRes_TruncErrorZero(iSpecies);

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

				/*--- Compute closest normal neighbor ---*/
				double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;

				double *Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
				cos_max = -1.0;
				for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
					jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
					scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
					for(iDim = 0; iDim < geometry->GetnDim(); iDim++) {
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
				// Point_Normal = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();

				/*--- Vector --> Velocity_corrected ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;


				for (iVar = 0; iVar < nVar; iVar ++)
					Res_Visc[iVar] = 0.0;

				/*--- Set the residual, truncation error and velocity value ---*/
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					node[iPoint]->SetVelocity_Old(Vector, iSpecies);
					node[iPoint]->SetVel_ResConv_Zero(iSpecies);
					node[iPoint]->SetVel_ResVisc_Zero(iSpecies);
					node[iPoint]->SetVelRes_TruncErrorZero(iSpecies);
				}

				Coord_i = geometry->node[iPoint]->GetCoord();
				Coord_j = geometry->node[Point_Normal]->GetCoord();

				dist_ij = 0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
				dist_ij = sqrt(dist_ij);

				U_i = node[iPoint]->GetSolution() ;
				U_j = node[Point_Normal]->GetSolution();

				CHeatTransfer[val_marker][iVertex] = 0.0;

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					Temperature = node[Point_Normal]->GetTemperature(iSpecies);
					Temperature_Gradient = (Twall - Temperature)/dist_ij;
					Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
					gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
					Gamma = config->GetSpecies_Gamma(iSpecies);
					Gamma_Minus_One = Gamma - 1.0;
					cp = (Gamma / Gamma_Minus_One) * gas_Constant;
					heat_flux_factor = cp * Viscosity/PRANDTL;

					Res_Visc[loc + nDim+1] = heat_flux_factor * Temperature_Gradient*Area;

					CHeatTransfer[val_marker][iVertex] -= heat_flux_factor * Temperature_Gradient;
				}
				node[iPoint]->SubtractRes_Visc(Res_Visc);  // SIGN CHECK


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
						rhoovisc = Density/(Viscosity+EPS); // rho over viscosity

						Jacobian_i[loc + nDim+1][loc + 0] = -factor*rhoovisc*theta*(phi_rho+phi*phi_p);
						Jacobian_i[loc + nDim+1][loc + nDim+1] = -factor*rhoovisc*(Gamma-1)*theta*phi_p;

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
					for(iDim = 0; iDim < geometry->GetnDim(); iDim++) {
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
				//	Point_Normal = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();

				/*--- Vector --> Velocity_corrected ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;


				for (iVar = 0; iVar < nVar; iVar ++)
					Res_Visc[iVar] = 0.0;

				/*--- Set the residual, truncation error and velocity value ---*/
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					node[iPoint]->SetVelocity_Old(Vector, iSpecies);
					node[iPoint]->SetVel_ResConv_Zero(iSpecies);
					node[iPoint]->SetVel_ResVisc_Zero(iSpecies);
					node[iPoint]->SetVelRes_TruncErrorZero(iSpecies);
				}

				Coord_i = geometry->node[iPoint]->GetCoord();
				Coord_j = geometry->node[Point_Normal]->GetCoord();

				dist_ij = 0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
				dist_ij = sqrt(dist_ij);

				U_i = node[iPoint]->GetSolution() ;
				U_j = node[Point_Normal]->GetSolution();

				CHeatTransfer[val_marker][iVertex] = 0.0;
				ion_density    = U_j[1*(nDim+2) + 0];
				atom_density   = U_j[0*(nDim+2) + 0];
				ionTemperature = node[Point_Normal]->GetTemperature(1);
				ion_flux       = - 0.25 * ion_density * catalyticity_coeff * sqrt( 8*Kb*ionTemperature/ ( PI_NUMBER* ionmass));
				electron_flux  = ion_flux*Mass[2]/Mass[1];
				atom_flux      = - ( ion_flux + electron_flux);

				Res_Visc[0*(nDim+2) + 0] = atom_flux     * Area;
				Res_Visc[1*(nDim+2) + 0] = ion_flux      * Area;
				Res_Visc[2*(nDim+2) + 0] = electron_flux * Area;

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					Temperature = node[Point_Normal]->GetTemperature(iSpecies);
					Temperature_Gradient = (Twall - Temperature)/dist_ij;
					Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
					gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
					Gamma = config->GetSpecies_Gamma(iSpecies);
					Gamma_Minus_One = Gamma - 1.0;
					cp = (Gamma / Gamma_Minus_One) * gas_Constant;
					heat_flux_factor = cp * Viscosity/PRANDTL;
					Res_Visc[loc + nDim+1] = heat_flux_factor * Temperature_Gradient*Area;
					CHeatTransfer[val_marker][iVertex] -= heat_flux_factor * Temperature_Gradient;
				}
				node[iPoint]->SubtractRes_Visc(Res_Visc);  // SIGN CHECK


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
						rhoovisc = Density/(Viscosity+EPS); // rho over viscosity

						Jacobian_i[loc + nDim+1][loc + 0] = -factor*rhoovisc*theta*(phi_rho+phi*phi_p);
						Jacobian_i[loc + nDim+1][loc + nDim+1] = -factor*rhoovisc*(Gamma-1)*theta*phi_p;

					}

					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				}
			}
		}
	}
}

void CPlasmaSolution::BC_Electrode(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

}

void CPlasmaSolution::BC_Dielectric(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

}


/*!
 * \method BC_Inlet
 * \brief Inlet Boundary Condition
 * \author A. Lonkar
 */

void CPlasmaSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
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


	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

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

			solver->GetPMatrix_inv(rho, velocity, c, kappa, invP_Matrix);
			solver->GetPMatrix(rho, velocity, c, kappa, P_Matrix);

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
			solver->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_update);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);

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

void CPlasmaSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

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
			solver->GetPMatrix_inv(rho, velocity, c, UnitaryNormal, invP_Matrix);
			solver->GetPMatrix(rho, velocity, c, UnitaryNormal, P_Matrix);

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

			solver->GetInviscidProjFlux(rho, velocity, pressure, enthalpy, UnitaryNormal, Residual);

			for	(iVar = 0; iVar < nVar; iVar++)
				Residual[iVar] = Residual[iVar]*Area;

			node[iPoint]->AddRes_Conv(Residual);
			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
				solver->SetConservative(U_domain, U_outlet);
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
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
void CPlasmaSolution::BC_Neumann(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint, Point_Normal, total_index;
	unsigned short iVar;
	double *U_domain, *U_interior;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	U_domain = new double[nVar];
	U_interior = new double[nVar];

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		Point_Normal = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			for (iVar = 0; iVar < nVar; iVar++) {
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);
				U_interior[iVar] = node[Point_Normal]->GetSolution(iVar);
			}
			node[iPoint]->SetSolution_Old(U_interior);
			node[iPoint]->SetSolution(U_interior);
			node[iPoint]->Set_ResConv_Zero();
			node[iPoint]->Set_ResSour_Zero();
			node[iPoint]->Set_ResVisc_Zero();
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
void CPlasmaSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, iSpecies, loc;
	unsigned short iVar, iDim;
	double *U_domain, *U_infty;

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	U_domain = new double[nVar];
	U_infty = new double[nVar];

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

				U_infty[loc + 0] = GetDensity_Inf(iSpecies);
				for (iDim = 0; iDim < nDim; iDim++)
					U_infty[loc+iDim+1] = GetDensity_Velocity_Inf(iDim, iSpecies);
				U_infty[loc+nDim+1] = GetDensity_Energy_Inf(iSpecies);
				if (iSpecies < nDiatomics)
					U_infty[loc+nDim+2] = GetDensity_Energy_vib_Inf(iSpecies);

				/*				U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iSpecies);
				U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iSpecies);
				U_infty[loc + 3] = GetDensity_Energy_Inf(iSpecies);
				if (nDim == 3) {
					U_infty[loc + 3] = GetDensity_Velocity_Inf(2,iSpecies);
					U_infty[loc + 4] = GetDensity_Energy_Inf(iSpecies);
				}*/
			}

			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_infty);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
	delete [] U_domain;
	delete [] U_infty;
}

/*!
 * \method BC_Sym_Plane
 * \brief Symmetry Boundary Condition
 * \author A. Lonkar
 */

void CPlasmaSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar, iSpecies, loc;
	double Pressure, Energy_el, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	double Area;
	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Compute the projected residual ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

			Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

				Residual[loc+0] = 0.0;
				Residual[loc+nDim+1] = 0.0;
				if ( iSpecies < nDiatomics )
					Residual[loc+nDim+2] = 0.0;
				Pressure = node[iPoint]->GetPressure(iSpecies);
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[loc + iDim+1] = Pressure*UnitaryNormal[iDim]*Area;
			}
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);


			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				double a2 = 0.0;
				double phi = 0.0;
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
          Energy_el = 0.0;
					a2 = config->GetSpecies_Gamma(iSpecies)-1.0;
					phi = a2*(0.5*node[iPoint]->GetVelocity2(iSpecies) - config->GetEnthalpy_Formation(iSpecies) - Energy_el);

					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

					for (iVar =  0; iVar <  nVar; iVar++) {
						Jacobian_i[loc + 0][iVar] = 0.0;
						Jacobian_i[loc + nDim+1][iVar] = 0.0;
					}
          if (iSpecies < nDiatomics) {
            for (iVar = 0; iVar < nVar; iVar++)
              Jacobian_i[loc+nDim+2][iVar] = 0.0;
          }
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
            if (iSpecies < nDiatomics)
              Jacobian_i[loc+iDim+1][loc+nDim+2] = a2*Normal[iDim];
					}
				}

				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}

		}
	}
}

void CPlasmaSolution::MPI_Send_Receive(CGeometry ***geometry, CSolution ****solution_container,
		CConfig **config, unsigned short iMGLevel, unsigned short iZone) {
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint;
	double *Conserv_Var, *Conserv_Undivided_Laplacian, 
	**PrimVar_Grad = NULL, **Conserv_Grad = NULL;

	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++)
		if (config[iZone]->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

			short SendRecv = config[iZone]->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry[iZone][iMGLevel]->nVertex[iMarker];
			bool viscous = config[iZone]->GetKind_ViscNumScheme_Plasma();

			/*--- Send information  ---*/
			if (SendRecv > 0) {

#ifndef NO_MPI

				unsigned long nBuffer_VectorVar = geometry[iZone][iMGLevel]->nVertex[iMarker]*nVar;
				unsigned long nBuffer_VectorSpe = geometry[iZone][iMGLevel]->nVertex[iMarker]*nSpecies;
				unsigned long nBuffer_Scalar = geometry[iZone][iMGLevel]->nVertex[iMarker];
				int send_to = SendRecv-1;

				double *Buffer_Send_U = new double[nBuffer_VectorVar];

				/*--- Upwind scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_UPWIND) {

					double *Buffer_Send_Ux = new double[nBuffer_VectorVar];
					double *Buffer_Send_Uy = new double[nBuffer_VectorVar];
					double *Buffer_Send_Uz = new double[nBuffer_VectorVar];

					/*--- Viscous part ---*/
					double *Buffer_Send_LaminarViscosity_S1  = NULL,
							*Buffer_Send_LaminarViscosity_S2 = NULL,
							*Buffer_Send_LaminarViscosity_S3 = NULL,
							*Buffer_Send_Vx = NULL,	*Buffer_Send_Vy = NULL, *Buffer_Send_Vz = NULL;
					if (viscous) {
						Buffer_Send_LaminarViscosity_S1 = new double[nBuffer_Scalar];
						Buffer_Send_LaminarViscosity_S2 = new double[nBuffer_Scalar];
						Buffer_Send_LaminarViscosity_S3 = new double[nBuffer_Scalar];
						Buffer_Send_Vx					= new double[nBuffer_VectorVar];
						Buffer_Send_Vy 					= new double[nBuffer_VectorVar];
						Buffer_Send_Vz 					= new double[nBuffer_VectorVar];
					}

					for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						Conserv_Var = node[iPoint]->GetSolution();
						if (iMGLevel == MESH_0) Conserv_Grad = node[iPoint]->GetGradient();

						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_U[iVar*nVertex+iVertex] = Conserv_Var[iVar];
							if (iMGLevel == MESH_0) {
								Buffer_Send_Ux[iVar*nVertex+iVertex] = Conserv_Grad[iVar][0];
								Buffer_Send_Uy[iVar*nVertex+iVertex] = Conserv_Grad[iVar][1];
								if (nDim == 3) Buffer_Send_Uz[iVar*nVertex+iVertex] = Conserv_Grad[iVar][2];
							}
						}
						if (viscous) {
							PrimVar_Grad = node[iPoint]->GetGradient_Primitive();
							for (iVar = 0; iVar < nVar; iVar++) {
								Buffer_Send_Vx[iVar*nVertex+iVertex] = PrimVar_Grad[iVar][0];
								Buffer_Send_Vy[iVar*nVertex+iVertex] = PrimVar_Grad[iVar][1];
								if (nDim == 3) Buffer_Send_Vz[iVar*nVertex+iVertex] = PrimVar_Grad[iVar][2];
							}
							Buffer_Send_LaminarViscosity_S1[iVertex] 	= node[iPoint]->GetLaminarViscosity(0);
							Buffer_Send_LaminarViscosity_S2[iVertex] 	= node[iPoint]->GetLaminarViscosity(1);
							Buffer_Send_LaminarViscosity_S3[iVertex] 	= node[iPoint]->GetLaminarViscosity(2);

						}
					}

					MPI::COMM_WORLD.Bsend(Buffer_Send_U,nBuffer_VectorVar,MPI::DOUBLE,send_to, 0);
					if (iMGLevel == MESH_0) {
						MPI::COMM_WORLD.Bsend(Buffer_Send_Ux,nBuffer_VectorVar,MPI::DOUBLE,send_to, 1);
						MPI::COMM_WORLD.Bsend(Buffer_Send_Uy,nBuffer_VectorVar,MPI::DOUBLE,send_to, 2);
						if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_Uz,nBuffer_VectorVar,MPI::DOUBLE,send_to, 3);
					}

					if (viscous) {
						MPI::COMM_WORLD.Bsend(Buffer_Send_LaminarViscosity_S1, nBuffer_Scalar, MPI::DOUBLE, send_to, 4);
						MPI::COMM_WORLD.Bsend(Buffer_Send_LaminarViscosity_S2, nBuffer_Scalar, MPI::DOUBLE, send_to, 5);
						MPI::COMM_WORLD.Bsend(Buffer_Send_LaminarViscosity_S3, nBuffer_Scalar, MPI::DOUBLE, send_to, 6);
						MPI::COMM_WORLD.Bsend(Buffer_Send_Vx, nBuffer_VectorVar, MPI::DOUBLE, send_to, 7);
						MPI::COMM_WORLD.Bsend(Buffer_Send_Vy, nBuffer_VectorVar, MPI::DOUBLE, send_to, 8);
						if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_Vz, nBuffer_VectorVar, MPI::DOUBLE, send_to, 9);
					}

					delete [] Buffer_Send_U;
					delete [] Buffer_Send_Ux;
					delete [] Buffer_Send_Uy;
					delete [] Buffer_Send_Uz;
					if (viscous) {
						delete [] Buffer_Send_Vx;
						delete [] Buffer_Send_Vy;
						delete [] Buffer_Send_Vz;
						delete [] Buffer_Send_LaminarViscosity_S1;
						delete [] Buffer_Send_LaminarViscosity_S2;
						delete [] Buffer_Send_LaminarViscosity_S3;
					}
				}

				/*--- Centered scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_CENTERED) {

					double *Buffer_Send_Undivided_Laplacian = NULL, *Buffer_Send_Sensor = NULL;
					Buffer_Send_Undivided_Laplacian = new double [nBuffer_VectorVar];
					Buffer_Send_Sensor = new double [nBuffer_Scalar];
					double *Buffer_Send_Lambda = new double[nBuffer_VectorSpe];
					unsigned short *Buffer_Send_Neighbor = new unsigned short [nBuffer_Scalar];

					/*--- Copy all the variables to the buffer vectors ---*/
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						Conserv_Var = node[iPoint]->GetSolution();
						Conserv_Undivided_Laplacian = node[iPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_U[iVar*nVertex+iVertex] = Conserv_Var[iVar];
							Buffer_Send_Undivided_Laplacian[iVar*nVertex+iVertex] = Conserv_Undivided_Laplacian[iVar];					
						}
						for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
							Buffer_Send_Lambda[iSpecies*nVertex+iVertex] = node[iPoint]->GetLambda(iSpecies);
						}
						Buffer_Send_Sensor[iVertex] = node[iPoint]->GetSensor();
						Buffer_Send_Neighbor[iVertex] = geometry[iZone][iMGLevel]->node[iPoint]->GetnPoint();
					}

					/*--- Send the buffer information ---*/
					MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_VectorVar, MPI::DOUBLE, send_to, 0);
					MPI::COMM_WORLD.Bsend(Buffer_Send_Undivided_Laplacian, nBuffer_VectorVar, MPI::DOUBLE, send_to, 1);
					MPI::COMM_WORLD.Bsend(Buffer_Send_Sensor, nBuffer_Scalar, MPI::DOUBLE, send_to, 2);
					MPI::COMM_WORLD.Bsend(Buffer_Send_Lambda, nBuffer_VectorSpe, MPI::DOUBLE, send_to, 3);
					MPI::COMM_WORLD.Bsend(Buffer_Send_Neighbor, nBuffer_Scalar, MPI::UNSIGNED_SHORT, send_to, 4);

					delete [] Buffer_Send_Undivided_Laplacian;
					delete [] Buffer_Send_Sensor;
					delete [] Buffer_Send_Lambda;
					delete [] Buffer_Send_Neighbor;
					delete [] Buffer_Send_U;

				}

#endif

			}

			/*--- Receive information  ---*/
			if (SendRecv < 0) {

				double rotMatrix[3][3], *angles;
				double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
				unsigned short iPeriodic_Index, iSpecies, loc;

				double *newSolution = new double[nVar];
				double **newGradient;
				newGradient = new double* [nVar];
				for (iVar = 0; iVar < nVar; iVar++)
					newGradient[iVar] = new double[3];

				unsigned long nBuffer_VectorVar = geometry[iZone][iMGLevel]->nVertex[iMarker]*nVar;
				unsigned long nBuffer_VectorSpe = geometry[iZone][iMGLevel]->nVertex[iMarker]*nSpecies;
				unsigned long nBuffer_Scalar = geometry[iZone][iMGLevel]->nVertex[iMarker];
				int receive_from = abs(SendRecv)-1;

				double *Buffer_Receive_U = new double[nBuffer_VectorVar];

				/*--- Upwind scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_UPWIND) {

					double *Buffer_Receive_Ux = new double[nBuffer_VectorVar];
					double *Buffer_Receive_Uy = new double[nBuffer_VectorVar];
					double *Buffer_Receive_Uz = new double[nBuffer_VectorVar];					

					/*--- Viscous part ---*/
					double  *Buffer_Receive_LaminarViscosity_S1 = NULL,	*Buffer_Receive_LaminarViscosity_S2 = NULL,
							*Buffer_Receive_LaminarViscosity_S3 = NULL,
							*Buffer_Receive_Vx = NULL, *Buffer_Receive_Vy = NULL, *Buffer_Receive_Vz = NULL;
					if (viscous) {
						Buffer_Receive_LaminarViscosity_S1 	= new double [nBuffer_Scalar];
						Buffer_Receive_LaminarViscosity_S2 	= new double [nBuffer_Scalar];
						Buffer_Receive_LaminarViscosity_S3 	= new double [nBuffer_Scalar];
						Buffer_Receive_Vx = new double [nBuffer_VectorVar];
						Buffer_Receive_Vy = new double [nBuffer_VectorVar];
						Buffer_Receive_Vz = new double [nBuffer_VectorVar];
					}


#ifdef NO_MPI
					/*--- Allow for periodic boundaries to use SEND_RECEIVE in serial.
					 Serial computations will only use the BC in receive mode, as
					 the proc will always be sending information to itself. ---*/

					/*--- Retrieve the donor boundary marker ---*/
					unsigned short donor_marker = 1;
					unsigned long donorPoint;
					for (unsigned short iMark = 0; iMark < config[iZone]->GetnMarker_All(); iMark++)
						if (config[iZone]->GetMarker_All_SendRecv(iMark)-1 == receive_from) donor_marker = iMark;

					/*--- Get the information from the donor point directly. This is a
					 serial computation with access to all nodes. Note that there is an
					 implicit ordering in the list. ---*/
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						donorPoint = geometry[iZone][iMGLevel]->vertex[donor_marker][iVertex]->GetNode();
						Conserv_Var = node[donorPoint]->GetSolution();
						Conserv_Grad = node[donorPoint]->GetGradient();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Receive_U[iVar*nVertex+iVertex] = Conserv_Var[iVar];
							Buffer_Receive_Ux[iVar*nVertex+iVertex] = Conserv_Grad[iVar][0];
							Buffer_Receive_Uy[iVar*nVertex+iVertex] = Conserv_Grad[iVar][1];
							if (nDim == 3) Buffer_Receive_Uz[iVar*nVertex+iVertex] = Conserv_Grad[iVar][2];
						}
						if (viscous) {
							PrimVar_Grad = node[donorPoint]->GetGradient_Primitive();
							for (iVar = 0; iVar < nVar; iVar++) {
								Buffer_Receive_Vx[iVar*nVertex+iVertex] = PrimVar_Grad[iVar][0];
								Buffer_Receive_Vy[iVar*nVertex+iVertex] = PrimVar_Grad[iVar][1];
								if (nDim == 3) Buffer_Receive_Vz[iVar*nVertex+iVertex] = PrimVar_Grad[iVar][2];
							}
							Buffer_Receive_LaminarViscosity_S1[iVertex] = node[donorPoint]->GetLaminarViscosity(0);
							Buffer_Receive_LaminarViscosity_S2[iVertex] = node[donorPoint]->GetLaminarViscosity(1);
							Buffer_Receive_LaminarViscosity_S3[iVertex] = node[donorPoint]->GetLaminarViscosity(2);

						}

					}

#else

					MPI::COMM_WORLD.Recv(Buffer_Receive_U,nBuffer_VectorVar,MPI::DOUBLE,receive_from, 0);
					if (iMGLevel == MESH_0) {
						MPI::COMM_WORLD.Recv(Buffer_Receive_Ux,nBuffer_VectorVar,MPI::DOUBLE,receive_from, 1);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Uy,nBuffer_VectorVar,MPI::DOUBLE,receive_from, 2);
						if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_Uz,nBuffer_VectorVar,MPI::DOUBLE,receive_from, 3);
					}

					if (viscous) {
						MPI::COMM_WORLD.Recv(Buffer_Receive_LaminarViscosity_S1, nBuffer_Scalar, MPI::DOUBLE, receive_from, 4);
						MPI::COMM_WORLD.Recv(Buffer_Receive_LaminarViscosity_S2, nBuffer_Scalar, MPI::DOUBLE, receive_from, 5);
						MPI::COMM_WORLD.Recv(Buffer_Receive_LaminarViscosity_S3, nBuffer_Scalar, MPI::DOUBLE, receive_from, 6);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Vx, nBuffer_VectorVar, MPI::DOUBLE, receive_from, 7);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Vy, nBuffer_VectorVar, MPI::DOUBLE, receive_from, 8);
						if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_Vz, nBuffer_VectorVar, MPI::DOUBLE, receive_from, 9);
					}
#endif

					/*--- Do the coordinate transformation ---*/
					for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {

						/*--- Find point and its type of transformation ---*/
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						iPeriodic_Index = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetRotation_Type();

						/*--- Retrieve the supplied periodic information. ---*/
						angles = config[iZone]->GetPeriodicRotation(iPeriodic_Index);

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

						/*--- Copy solution before performing transformation. ---*/
						for (iVar = 0; iVar < nVar; iVar++) newSolution[iVar] = Buffer_Receive_U[iVar*nVertex+iVertex];

						/*--- Need to rotate the momentum components. ---*/
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
							if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

							if (nDim == 2) {
								newSolution[loc+1] = rotMatrix[0][0]*Buffer_Receive_U[(loc+1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[(loc+2)*nVertex+iVertex];
								newSolution[loc+2] = rotMatrix[1][0]*Buffer_Receive_U[(loc+1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[(loc+2)*nVertex+iVertex];
							}
							else {
								newSolution[loc+1] = rotMatrix[0][0]*Buffer_Receive_U[(loc+1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[(loc+2)*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_U[(loc+3)*nVertex+iVertex];
								newSolution[loc+2] = rotMatrix[1][0]*Buffer_Receive_U[(loc+1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[(loc+2)*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_U[(loc+3)*nVertex+iVertex];
								newSolution[loc+3] = rotMatrix[2][0]*Buffer_Receive_U[(loc+1)*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_U[(loc+2)*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_U[(loc+3)*nVertex+iVertex];
							}
						}

						/*--- Copy transformed conserved variables back into buffer. ---*/
						for (iVar = 0; iVar < nVar; iVar++) Buffer_Receive_U[iVar*nVertex+iVertex] = newSolution[iVar];

						/*--- Also transform the gradient for upwinding if this is the fine mesh ---*/
						for (iVar = 0; iVar < nVar; iVar++) {
							newGradient[iVar][0] = Buffer_Receive_Ux[iVar*nVertex+iVertex];
							newGradient[iVar][1] = Buffer_Receive_Uy[iVar*nVertex+iVertex];
							if (nDim == 3) newGradient[iVar][2] = Buffer_Receive_Uz[iVar*nVertex+iVertex];
						}

						/*--- Need to rotate the gradients for all conserved variables. ---*/
						for (iVar = 0; iVar < nVar; iVar++) {
							if (nDim == 2) {
								newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Ux[iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_Uy[iVar*nVertex+iVertex];
								newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Ux[iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_Uy[iVar*nVertex+iVertex];
							}
							else {
								newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Ux[iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_Uy[iVar*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_Uz[iVar*nVertex+iVertex];
								newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Ux[iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_Uy[iVar*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_Uz[iVar*nVertex+iVertex];
								newGradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Ux[iVar*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_Uy[iVar*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_Uz[iVar*nVertex+iVertex];
							}
						}

						/*--- Copy transformed gradients back into buffer. ---*/
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Receive_Ux[iVar*nVertex+iVertex] = newGradient[iVar][0];
							Buffer_Receive_Uy[iVar*nVertex+iVertex] = newGradient[iVar][1];
							if (nDim == 3) Buffer_Receive_Uz[iVar*nVertex+iVertex] = newGradient[iVar][2];
						}

						if (viscous) {
							for (iVar = 0; iVar < nVar; iVar++) {
								newGradient[iVar][0] = Buffer_Receive_Vx[iVar*nVertex+iVertex];
								newGradient[iVar][1] = Buffer_Receive_Vy[iVar*nVertex+iVertex];
								if (nDim == 3) newGradient[iVar][2] = Buffer_Receive_Vz[iVar*nVertex+iVertex];
							}

							/*--- Need to rotate the gradients for all conserved variables. ---*/
							for (iVar = 0; iVar < nVar; iVar++) {
								if (nDim == 2) {
									newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Vx[iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_Vy[iVar*nVertex+iVertex];
									newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Vx[iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_Vy[iVar*nVertex+iVertex];
								}
								else {
									newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Vx[iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_Vy[iVar*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_Vz[iVar*nVertex+iVertex];
									newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Vx[iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_Vy[iVar*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_Vz[iVar*nVertex+iVertex];
									newGradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Vx[iVar*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_Vy[iVar*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_Vz[iVar*nVertex+iVertex];
								}
							}

							/*--- Copy transformed gradients back into buffer. ---*/
							for (iVar = 0; iVar < nVar; iVar++) {
								Buffer_Receive_Vx[iVar*nVertex+iVertex] = newGradient[iVar][0];
								Buffer_Receive_Vy[iVar*nVertex+iVertex] = newGradient[iVar][1];
								if (nDim == 3) Buffer_Receive_Vz[iVar*nVertex+iVertex] = newGradient[iVar][2];
							}
						}

						/*--- Store the received information ---*/
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
							node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Ux[iVar*nVertex+iVertex]);
							node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Uy[iVar*nVertex+iVertex]);
							if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Uz[iVar*nVertex+iVertex]);
						}

						if (viscous) {
							for (iVar = 0; iVar < nVar; iVar++) {
								node[iPoint]->SetGradient_Primitive(iVar, 0, Buffer_Receive_Vx[iVar*nVertex+iVertex]);
								node[iPoint]->SetGradient_Primitive(iVar, 1, Buffer_Receive_Vy[iVar*nVertex+iVertex]);
								if (nDim == 3) node[iPoint]->SetGradient_Primitive(iVar, 2, Buffer_Receive_Vz[iVar*nVertex+iVertex]);
							}
							node[iPoint]->SetLaminarViscosity(Buffer_Receive_LaminarViscosity_S1[iVertex], 0);
							node[iPoint]->SetLaminarViscosity(Buffer_Receive_LaminarViscosity_S2[iVertex], 1);
							node[iPoint]->SetLaminarViscosity(Buffer_Receive_LaminarViscosity_S3[iVertex], 2);
						}

					}

					delete [] Buffer_Receive_U;
					delete [] Buffer_Receive_Ux;
					delete [] Buffer_Receive_Uy;
					delete [] Buffer_Receive_Uz;
					if (viscous) {
						delete [] Buffer_Receive_LaminarViscosity_S1;
						delete [] Buffer_Receive_LaminarViscosity_S2;
						delete [] Buffer_Receive_LaminarViscosity_S3;
						delete [] Buffer_Receive_Vx;
						delete [] Buffer_Receive_Vy;
						delete [] Buffer_Receive_Vz;
					}

					delete [] newSolution;
					for (iVar = 0; iVar < nVar; iVar++)
						delete [] newGradient[iVar];
					delete [] newGradient;
				}

				/*--- Centered scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_CENTERED) {

					double *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Receive_Sensor = NULL;
					Buffer_Receive_Undivided_Laplacian = new double [nBuffer_VectorVar];
					Buffer_Receive_Sensor = new double [nBuffer_Scalar];
					double *Buffer_Receive_Lambda = new double [nBuffer_VectorSpe];
					unsigned short *Buffer_Receive_Neighbor = new unsigned short [nBuffer_Scalar];


#ifdef NO_MPI
					/*--- Allow for periodic boundaries to use SEND_RECEIVE in serial.
					 Serial computations will only use the BC in receive mode, as
					 the proc will always be sending information to itself. ---*/

					/*--- Retrieve the donor information from the matching marker ---*/
					unsigned short donor_marker = 1;
					unsigned long donorPoint;
					for (unsigned short iMark = 0; iMark < config[iZone]->GetnMarker_All(); iMark++)
						if (config[iZone]->GetMarker_All_SendRecv(iMark)-1 == receive_from) donor_marker = iMark;				

					/*--- Get the information from the donor directly. This is a serial
					 computation with access to all nodes. Note that there is an
					 implicit ordering in the list. ---*/
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						donorPoint = geometry[iZone][iMGLevel]->vertex[donor_marker][iVertex]->GetNode();
						Conserv_Var = node[donorPoint]->GetSolution();
						Conserv_Undivided_Laplacian = node[donorPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Receive_U[iVar*nVertex+iVertex] = Conserv_Var[iVar];
							Buffer_Receive_Undivided_Laplacian[iVar*nVertex+iVertex] = Conserv_Undivided_Laplacian[iVar];					
						}
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
							Buffer_Receive_Lambda[iSpecies*nVertex+iVertex] = node[donorPoint]->GetLambda();
						}
						Buffer_Receive_Sensor[iVertex] = node[donorPoint]->GetSensor();
						Buffer_Receive_Neighbor[iVertex] = geometry[iZone][iMGLevel]->node[donorPoint]->GetnPoint();
					}
#else
					/*--- Receive the information ---*/
					MPI::COMM_WORLD.Recv(Buffer_Receive_U,nBuffer_VectorVar,MPI::DOUBLE,receive_from, 0);
					MPI::COMM_WORLD.Recv(Buffer_Receive_Undivided_Laplacian,nBuffer_VectorVar,MPI::DOUBLE,receive_from, 1);
					MPI::COMM_WORLD.Recv(Buffer_Receive_Sensor,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
					MPI::COMM_WORLD.Recv(Buffer_Receive_Lambda,nBuffer_VectorSpe,MPI::DOUBLE,receive_from, 3);
					MPI::COMM_WORLD.Recv(Buffer_Receive_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,receive_from, 4);
#endif

					/*--- Do the coordinate transformation ---*/
					for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {

						/*--- Find point and its type of transformation ---*/
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						iPeriodic_Index = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetRotation_Type();

						/*--- Retrieve the supplied periodic information. ---*/
						angles = config[iZone]->GetPeriodicRotation(iPeriodic_Index);

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

						/*--- Copy solution before performing transformation. ---*/
						for (iVar = 0; iVar < nVar; iVar++) newSolution[iVar] = Buffer_Receive_U[iVar*nVertex+iVertex];

						/*--- Need to rotate the momentum components. ---*/
						if (nDim == 2) {
							newSolution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[2*nVertex+iVertex];
							newSolution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[2*nVertex+iVertex];
						} 
						else {
							newSolution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[2*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_U[3*nVertex+iVertex];
							newSolution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[2*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_U[3*nVertex+iVertex];
							newSolution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_U[2*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_U[3*nVertex+iVertex];
						}

						/*--- Copy transformed conserved variables back into buffer. ---*/
						for (iVar = 0; iVar < nVar; iVar++) Buffer_Receive_U[iVar*nVertex+iVertex] = newSolution[iVar];

						/*--- Store the received information ---*/
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
							node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVar*nVertex+iVertex]);
						}

						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
							node[iPoint]->SetLambda(Buffer_Receive_Lambda[iSpecies*nVertex+iVertex], iSpecies);
						}

						node[iPoint]->SetSensor(Buffer_Receive_Sensor[iVertex]);

						geometry[iZone][iMGLevel]->node[iPoint]->SetnNeighbor(Buffer_Receive_Neighbor[iVertex]);

					}

					delete [] Buffer_Receive_Undivided_Laplacian;
					delete [] Buffer_Receive_Sensor;
					delete [] Buffer_Receive_Lambda;
					delete [] Buffer_Receive_Neighbor;
					delete [] Buffer_Receive_U;

					delete [] newSolution;
					for (iVar = 0; iVar < nVar; iVar++)
						delete [] newGradient[iVar];
					delete [] newGradient;

				}	

			}
		}

}
