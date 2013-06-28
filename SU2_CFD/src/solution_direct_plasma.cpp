/*!
 * \file solution_direct_plasma.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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
	double Vel2 = 0.0, Ru;
	string mesh_filename, text_line;
	ifstream restart_file;

	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	if (config->GetKind_GasModel() == ARGON) {
		Gamma = config->GetGamma();
		Gamma_Minus_One = Gamma - 1.0;

		/*--- Define geometry constans in the solver structure ---*/
		nDim = geometry->GetnDim();
		nFluids = config->GetnFluids() ;
		nSpecies = config->GetnSpecies();
		nMonatomics = nSpecies;
		nDiatomics = 0;
		nVar =  nSpecies + nFluids*nDim + nFluids;	// Continuity + nDim*Momentum + Energy;
		node = new CVariable*[geometry->GetnPoint()];

		double M1, M2, M3;
		double R1, R2, R3;


		/*--- Define  some auxiliar vector related with the residual ---*/
		Residual = new double[nVar];	Residual_Max = new double[nVar];
		Residual_i = new double[nVar];	Residual_j = new double[nVar];
		Res_Conv = new double[nVar];	Res_Visc = new double[nVar];	Res_Sour = new double[nVar];

		/*--- Define some auxiliar vector related with the solution ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];

		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector = new double[nDim];
		Vector_i = new double[nDim]; Vector_j = new double[nDim];

		/*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
		if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) {
			p1_Und_Lapl = new double [geometry->GetnPoint()]; p2_Und_Lapl = new double [geometry->GetnPoint()]; }

		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT) {
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar]; Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar]; Jacobian_j[iVar] = new double [nVar]; }
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Initialize_Jacobian_Structure(geometry, config);
			xsol = new double [geometry->GetnPoint()*nVar];
			rhs = new double [geometry->GetnPoint()*nVar];
		}

		/*--- Computation of gradients by least squares ---*/
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
				(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
			/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
			Smatrix = new double* [nDim];
			for (iDim = 0; iDim < nDim; iDim++)
				Smatrix[iDim] = new double [nDim];
			/*--- c vector := transpose(WA)*(Wb) ---*/
			cvector = new double* [nVar+nSpecies];
			for (iVar = 0; iVar < nVar+nSpecies; iVar++)
				cvector[iVar] = new double [nDim];
		}

		CSkinFriction = new double* [config->GetnMarker_All()];
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];

		CHeatTransfer = new double* [config->GetnMarker_All()];
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			CHeatTransfer[iMarker] = new double [geometry->nVertex[iMarker]];

		Prandtl_Lam   = config->GetPrandtl_Lam();

		M1 = AVOGAD_CONSTANT*config->GetParticle_Mass(0);
		M2 = AVOGAD_CONSTANT*config->GetParticle_Mass(1);
		M3 = AVOGAD_CONSTANT*config->GetParticle_Mass(2);

		Ru = UNIVERSAL_GAS_CONSTANT; R1 = Ru/M1; R2 = Ru/M2; R3 = Ru/M3;

		double *Cv =  new double [nSpecies];

		Cv[0]  = R1/(Gamma-1.0);
		Cv[1]  = R2/(Gamma-1.0);
		Cv[2]  = R3/(Gamma-1.0);

		/*--- Flow infinity initialization stuff ---*/
		Density_Inf = new double [nSpecies];
		Pressure_Inf = new double [nSpecies];
		Velocity_Inf = new double*[nSpecies];
		Mach_Inf = 	new double [nSpecies];
		Energy_Inf = new double [nSpecies];
		Temperature_Inf = new double [nSpecies];


		Density_Inlet = new double [nSpecies];
		Velocity_Inlet  = new double*[nSpecies];
		Mach_Inlet 		= new double [nSpecies];
		Pressure_Inlet  = new double [nSpecies];
		Energy_Inlet = new double [nSpecies];
		Density_Outlet = new double [nSpecies];
		Pressure_Outlet = new double [nSpecies];
		Energy_Outlet = new double [nSpecies];
		Mach_Outlet = new double [nSpecies];
		Velocity_Outlet  = new double*[nSpecies];

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Velocity_Inf[iSpecies] = new double [nDim];
			Velocity_Inlet[iSpecies] = new double [nDim];
			Velocity_Outlet[iSpecies] = new double [nDim];
		}

#ifdef OneDArgonSimulation

		Density_Inf[0] = 2.1331E-1;
		Density_Inf[1] = 1E-3*Density_Inf[0]*M2/M1;
		Density_Inf[2] = 1E-3*Density_Inf[0]*M3/M1;

		Pressure_Inf[0] = .1335E5;
		Pressure_Inf[1] = .1335E2;
		Pressure_Inf[2] = .1335E2;

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Velocity_Inf[iSpecies][0] = 4800.0;
			for (iDim = 1; iDim < nDim; iDim ++)
				Velocity_Inf[iSpecies][iDim] = 0.0;
		}

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));

		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];
			}
			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*Gamma_Minus_One)+0.5*Vel2;
		}


		/*--- Flow at the inlet initialization stuff ---*/

		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
			Density_Inlet[iSpecies] = Density_Inf[iSpecies];

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Mach_Inlet[iSpecies]     = Mach_Inf[iSpecies];
			Pressure_Inlet[iSpecies] = Pressure_Inf[iSpecies];
			Energy_Inlet[iSpecies]   = Energy_Inf[iSpecies];
			for (iDim = 0; iDim < nDim; iDim ++)
				Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
		}

		/*--- Flow at the Outlet initialization stuff ---*/

		Density_Outlet[0] = 0.8745E0;
		Density_Outlet[1] = 1E-3*Density_Outlet[0]*M2/M1;
		Density_Outlet[2] = Density_Outlet[1]*M3/M2;

		Pressure_Outlet[0] = 0.6167E7;
		Pressure_Outlet[1] = 0.6149E4;
		Pressure_Outlet[2] = 0.5461E2;

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Velocity_Outlet[iSpecies][0] = 188.0;
			for (iDim = 1; iDim < nDim; iDim ++)
				Velocity_Outlet[iSpecies][iDim] = 0.0;
			Vel2 = 0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Vel2 += Velocity_Outlet[iSpecies][iDim]*Velocity_Outlet[iSpecies][iDim];
			}
			Energy_Outlet[iSpecies] = Pressure_Outlet[iSpecies]/(Density_Outlet[iSpecies]*Gamma_Minus_One)+0.5*Vel2;
			Mach_Outlet[iSpecies] = Velocity_Outlet[iSpecies][0]/sqrt(Gamma*Pressure_Outlet[iSpecies]/Density_Outlet[iSpecies]);

		}

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			if(geometry->node[iPoint]->GetCoord(0) <= 0.0)
				node[iPoint] = new CPlasmaVariable(Density_Inlet, Velocity_Inlet, Energy_Inlet,nDim, nVar, nSpecies, nFluids, config);
			else
				node[iPoint] = new CPlasmaVariable(Density_Outlet, Velocity_Outlet, Energy_Outlet, nDim, nVar, nSpecies, nFluids, config);
		}
# endif


#ifdef Debug_Cylinder_From_Infinity

		double usqr;
		double ionisation = 1;
		Temperature_Inf[0] = 273.15;
		Temperature_Inf[1] = 273.15;
		Temperature_Inf[2] = 273.15;


		Density_Inf[0] = 6.92252E-6;
		Density_Inf[1] = Density_Inf[0]* ionisation *M2/M1;
		Density_Inf[2] = Density_Inf[0]* ionisation *M3/M1;

		Pressure_Inf[0] = Density_Inf[0]*R1*Temperature_Inf[0];
		Pressure_Inf[1] = Density_Inf[1]*R2*Temperature_Inf[1];
		Pressure_Inf[2] = Density_Inf[2]*R3*Temperature_Inf[2];

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
			Velocity_Inf[iSpecies][0] = 99.5369;
			for (iDim = 1; iDim < nDim; iDim ++)
				Velocity_Inf[iSpecies][iDim] = 0.0;
			usqr = Velocity_Inf[iSpecies][0] *Velocity_Inf[iSpecies][0] ;

		}


		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Energy_Inf[iSpecies] = Cv[iSpecies]*Temperature_Inf[iSpecies] + 0.5*usqr;
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
		}

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Density_Outlet[iSpecies] 	= Density_Inf[iSpecies];
			Density_Inlet[iSpecies] 	= Density_Inf[iSpecies];
			Pressure_Outlet[iSpecies]   = Pressure_Inf[iSpecies];
			Pressure_Inlet[iSpecies]   = Pressure_Inf[iSpecies];
			Mach_Outlet[iSpecies] 		= Mach_Inf[iSpecies];
			Mach_Inlet[iSpecies] 		= Mach_Inf[iSpecies];
			Energy_Outlet[iSpecies] 	= Energy_Inf[iSpecies];
			Energy_Inlet[iSpecies] 	= Energy_Inf[iSpecies];

			for (iDim = 0; iDim < nDim; iDim ++) {
				Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
			}
		}

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			cout << " iSpecies = " << iSpecies << endl;
			cout << " Temperature_inf = " << Temperature_Inf[iSpecies] << endl;
			cout << " Pressure_Inf = " << Pressure_Inf[iSpecies] << " Pressure_Outlet = " << Pressure_Outlet[iSpecies] << endl;
			cout << " Density_Inf = " << Density_Inf[iSpecies] << " Density_Outlet = " << Density_Outlet[iSpecies] << endl;
			cout << " Mach_Inf = " << Mach_Inf[iSpecies] << " Mach_Outlet = " << Mach_Outlet[iSpecies] << endl;
			cout << " *********" << endl;
		}

		cout << " Gas_Constant = " << R1 << "  R2 = "  << R2 << " R3 =" << R3 << endl;

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, nSpecies, nFluids, config);

#endif



		//#ifdef Bullet_From_Infinity
		if (!restart) {
			double ionisation = 0.01;
			Temperature_Inf[0] = 810.0;
			Temperature_Inf[1] = 810.0;
			Temperature_Inf[2] = 3900.0;
			double uscalar = 2437.0;

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
				Density_Outlet[iSpecies] 	= Density_Inf[iSpecies];
				Density_Inlet[iSpecies] 	= Density_Inf[iSpecies];
				Pressure_Outlet[iSpecies]   = Pressure_Inf[iSpecies];
				Pressure_Inlet[iSpecies]   = Pressure_Inf[iSpecies];
				Mach_Outlet[iSpecies] 		= Mach_Inf[iSpecies];
				Mach_Inlet[iSpecies] 		= Mach_Inf[iSpecies];
				Energy_Outlet[iSpecies] 	= Energy_Inf[iSpecies];
				Energy_Inlet[iSpecies] 	= Energy_Inf[iSpecies];

				for (iDim = 0; iDim < nDim; iDim ++) {
					Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
					Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				}
			}

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				cout << " iSpecies = " << iSpecies << endl;
				cout << " Temperature_inf = " << Temperature_Inf[iSpecies] << endl;
				cout << " Pressure_Inf = " << Pressure_Inf[iSpecies] << " Pressure_Outlet = " << Pressure_Outlet[iSpecies] << endl;
				cout << " Density_Inf = " << Density_Inf[iSpecies] << " Density_Outlet = " << Density_Outlet[iSpecies] << endl;
				cout << " Mach_Inf = " << Mach_Inf[iSpecies] << " Mach_Outlet = " << Mach_Outlet[iSpecies] << endl;
				cout << " *********" << endl;
			}

			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, nSpecies, nFluids,nMonatomics, nDiatomics, config);
		}
		//#endif


#ifdef Bullet_EulerRestart
		if (restart) {
			string mesh_filename = config->GetSolution_FlowFileName();
			ifstream restart_file;

			char *cstr; cstr = new char [mesh_filename.size()+1];
			strcpy (cstr, mesh_filename.c_str());
			restart_file.open(cstr, ios::in);
			if (restart_file.fail()) {
				cout << "There is no flow restart file!!" << endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
				exit(1);
			}

			unsigned long index;
			string text_line;
			unsigned short nVar_Argon = nDim + 2;
			double ionisation = 0.01;
			double *u = new double[nDim];
			double **Velocity = new double*[nSpecies];
			double *Temperature = new double [nSpecies];
			double *Density = new double[nSpecies];
			double *Energy = new double [nSpecies];

			Temperature_Inf[0] = 810.0;
			Temperature_Inf[1] = 810.0;
			Temperature_Inf[2] = 3900.0;

			Density_Inf[0] = 1.662E-4;
			Density_Inf[1] = Density_Inf[0]* ionisation *M2/M1;
			Density_Inf[2] = Density_Inf[0]* ionisation *M3/M1;

			Pressure_Inf[0] = Density_Inf[0]*R1*Temperature_Inf[0];
			Pressure_Inf[1] = Density_Inf[1]*R2*Temperature_Inf[1];
			Pressure_Inf[2] = Density_Inf[2]*R3*Temperature_Inf[2];

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
				Velocity_Inf[iSpecies] = new double [nDim];
				Velocity_Inf[iSpecies][0] = 2437.0;
				for (iDim = 1; iDim < nDim; iDim ++)
					Velocity_Inf[iSpecies][iDim] = 0.0;
			}
			double usqr = 2437.0*2437.0;

			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				Energy_Inf[iSpecies] = Cv[iSpecies]*Temperature_Inf[iSpecies] + 0.5*usqr;
				Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			}

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				Density_Outlet[iSpecies] 	= Density_Inf[iSpecies];
				Density_Inlet[iSpecies] 	= Density_Inf[iSpecies];
				Pressure_Outlet[iSpecies]   = Pressure_Inf[iSpecies];
				Pressure_Inlet[iSpecies]    = Pressure_Inf[iSpecies];
				Mach_Outlet[iSpecies] 		= Mach_Inf[iSpecies];
				Mach_Inlet[iSpecies] 		= Mach_Inf[iSpecies];
				Energy_Outlet[iSpecies] 	= Energy_Inf[iSpecies];
				Energy_Inlet[iSpecies] 		= Energy_Inf[iSpecies];

				for (iDim = 0; iDim < nDim; iDim ++) {
					Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
					Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				}
			}

			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				getline(restart_file,text_line);
				istringstream point_line(text_line);
				point_line >> index;
				for (iVar = 0; iVar < nVar_Argon; iVar ++) {
					point_line >> Solution[iVar];
				}
				usqr = 0.0;
				for (iDim = 0; iDim < nDim; iDim ++) {
					u[iDim] = Solution[1+iDim]/Solution[0];
					usqr = usqr + u[iDim]*u[iDim];
				}

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					Velocity[iSpecies] = new double [nDim];
					for (iDim = 0; iDim < nDim; iDim ++)
						Velocity[iSpecies][iDim] = u[iDim];
				}

				Temperature[0] = (Solution[nDim+1]/Solution[0] - 0.5*usqr)/Cv[0];
				Temperature[1] = Temperature[0];
				Temperature[2] = Temperature_Inf[2];

				Density[0] = Solution[0];
				Density[1] = Solution[0] * ionisation * M2/M1;
				Density[2] = Solution[0] * ionisation * M3/M1;

				for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
					Energy[iSpecies] = Cv[iSpecies]*Temperature[iSpecies] + 0.5*usqr;

				node[iPoint] = new CPlasmaVariable(Density, Velocity, Energy, nDim, nVar, nSpecies, nFluids, nMonatomics, nDiatomics, config);
			}
			restart_file.close();
		}
#endif


		//#ifdef Bullet_PlasmaRestart
		if (restart) {
			string mesh_filename = config->GetSolution_FlowFileName();
			ifstream restart_file;

			char *cstr; cstr = new char [mesh_filename.size()+1];
			strcpy (cstr, mesh_filename.c_str());
			restart_file.open(cstr, ios::in);
			if (restart_file.fail()) {
				cout << "There is no flow restart file!!" << endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
				exit(1);
			}

			unsigned long index;
			string text_line;
//			unsigned short nVar_Argon = nDim + 2;
			double ionisation = 0.01;
//			double *u = new double[nDim];
//			double **Velocity = new double*[nSpecies];
//			double *Temperature = new double [nSpecies];
//			double *Density = new double[nSpecies];
//			double *Energy = new double [nSpecies];
			double uscalar = 2437.0;


			Temperature_Inf[0] = 810.0;
			Temperature_Inf[1] = 810.0;
			Temperature_Inf[2] = 3900.0;

			Density_Inf[0] = 1.662E-4;
			Density_Inf[1] = Density_Inf[0]* ionisation *M2/M1;
			Density_Inf[2] = Density_Inf[0]* ionisation *M3/M1;

			Pressure_Inf[0] = Density_Inf[0]*R1*Temperature_Inf[0];
			Pressure_Inf[1] = Density_Inf[1]*R2*Temperature_Inf[1];
			Pressure_Inf[2] = Density_Inf[2]*R3*Temperature_Inf[2];

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
				Velocity_Inf[iSpecies] = new double [nDim];
				Velocity_Inf[iSpecies][0] = uscalar;
				for (iDim = 1; iDim < nDim; iDim ++)
					Velocity_Inf[iSpecies][iDim] = 0.0;
			}
			double usqr = uscalar*uscalar;

			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				Energy_Inf[iSpecies] = Cv[iSpecies]*Temperature_Inf[iSpecies] + 0.5*usqr;
				Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			}

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				Density_Outlet[iSpecies] 	= Density_Inf[iSpecies];
				Density_Inlet[iSpecies] 	= Density_Inf[iSpecies];
				Pressure_Outlet[iSpecies]   = Pressure_Inf[iSpecies];
				Pressure_Inlet[iSpecies]    = Pressure_Inf[iSpecies];
				Mach_Outlet[iSpecies] 		= Mach_Inf[iSpecies];
				Mach_Inlet[iSpecies] 		= Mach_Inf[iSpecies];
				Energy_Outlet[iSpecies] 	= Energy_Inf[iSpecies];
				Energy_Inlet[iSpecies] 		= Energy_Inf[iSpecies];

				for (iDim = 0; iDim < nDim; iDim ++) {
					Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
					Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				}
			}
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				getline(restart_file,text_line);
				istringstream point_line(text_line);
				point_line >> index;
				for (iVar = 0; iVar < nVar; iVar ++)
					point_line >> Solution[iVar];
				node[iPoint] = new CPlasmaVariable(Solution, nDim, nVar, nSpecies, nFluids, nMonatomics, nDiatomics, config);

			}

			restart_file.close();
		}
		//#endif
	}
	if (config->GetKind_GasModel() == AIR7) {

		/*--- Define geometry constants in the solver structure ---*/
		nDim = geometry->GetnDim();
		nMonatomics = 7;
		nDiatomics = 0;
		nSpecies = nMonatomics + nDiatomics;
		nVar = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);
		node = new CVariable*[geometry->GetnPoint()];

		double GammaDiatomic = config->GetGammaDiatomic();
		double GammaMonatomic = config->GetGammaMonatomic();

		/*--- Define auxiliary vectors for residuals at nodes i & j ---*/
		Residual = new double[nVar];	Residual_Max = new double[nVar];
		Residual_Chemistry = new double[nVar]; Residual_MomentumExch = new double[nVar];
		Residual_ElecForce = new double[nVar]; Residual_EnergyExch = new double[nVar];
		Residual_i = new double[nVar];	Residual_j = new double[nVar];
		Res_Conv = new double[nVar];	Res_Visc = new double[nVar];	Res_Sour = new double[nVar];

		/*--- Define auxiliary vectors for the solution at nodes i & j ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];

		/*--- Define auxiliary vectors for the geometry ---*/
		Vector = new double[nDim];
		Vector_i = new double[nDim]; Vector_j = new double[nDim];

		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT) {
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar]; Jacobian_j = new double* [nVar];
			Jacobian_Chemistry = new double* [nVar];
			Jacobian_ElecForce = new double* [nVar];
			Jacobian_MomentumExch = new double* [nVar];
			Jacobian_EnergyExch = new double* [nVar];

			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar]; Jacobian_j[iVar] = new double [nVar]; 
				Jacobian_Chemistry[iVar] = new double [nVar];
				Jacobian_ElecForce[iVar] = new double [nVar];
				Jacobian_MomentumExch[iVar] = new double [nVar];
				Jacobian_EnergyExch[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Initialize_Jacobian_Structure(geometry, config);
			xsol = new double [geometry->GetnPoint()*nVar];
			rhs = new double [geometry->GetnPoint()*nVar];
		}

		/*--- Computation of gradients by least squares ---*/
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
				(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
			/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
			Smatrix = new double* [nDim];
			for (iDim = 0; iDim < nDim; iDim++)
				Smatrix[iDim] = new double [nDim];
			/*--- c vector := transpose(WA)*(Wb) ---*/
			cvector = new double* [nVar+nSpecies];
			for (iVar = 0; iVar < nVar+nSpecies; iVar++)
				cvector[iVar] = new double [nDim];
		}

		/*--- Flow infinity array initialization ---*/
		Density_Inf = new double [nSpecies];
		Pressure_Inf = new double [nSpecies];
		Velocity_Inf = new double*[nSpecies];
		Mach_Inf = 	new double [nSpecies];
		Energy_Inf = new double [nSpecies];
		Energy_vib_Inf = new double [nDiatomics];
		Enthalpy_Formation = new double [nSpecies];
		Molar_Mass = new double[nSpecies];

		Density_Inlet = new double [nSpecies];
		Velocity_Inlet  = new double*[nSpecies];
		Mach_Inlet 		= new double [nSpecies];
		Pressure_Inlet  = new double [nSpecies];
		Energy_Inlet = new double [nSpecies];
		Density_Outlet = new double [nSpecies];
		Pressure_Outlet = new double [nSpecies];
		Energy_Outlet = new double [nSpecies];
		Mach_Outlet = new double [nSpecies];
		Velocity_Outlet  = new double*[nSpecies];		
		Gas_Composition = new double[nSpecies];
		Mach2Vel_FreeStream = new double[nSpecies];

		/*--- Get Mean Flow Properties from the configuration file ---*/
		double Pressure_Inf_mean = config->GetPressure_FreeStream();
		double Temperature_Inf_mean = config->GetTemperature_FreeStream();
		double Gas_Constant_mean = config->GetGas_Constant();
		double Mach2Vel_FreeStream_mean = sqrt(config->GetGamma()*Gas_Constant_mean*Temperature_Inf_mean);
		double *Velocity_Inf_mean;
		Velocity_Inf_mean = new double[nDim];
		for (iDim = 0; iDim < nDim; iDim++) {
			if (iDim == 0)
				Velocity_Inf_mean[iDim] = config->GetMach_FreeStreamND()*Mach2Vel_FreeStream_mean;
			else 
				Velocity_Inf_mean[iDim] = 0.0;
		}

		/*--- Initialize species quanitites from mean flow properties ---*/
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			Enthalpy_Formation[iSpecies] = config->GetEnthalpy_Formation(iSpecies);
			Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies);
			Gas_Composition[iSpecies] = config->GetInitial_Gas_Composition(iSpecies);
		}
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			Velocity_Inf[iSpecies] = new double [nDim];
			Velocity_Inlet[iSpecies] = new double [nDim];
			Velocity_Outlet[iSpecies] = new double [nDim];

			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity_Inf[iSpecies][iDim] = Velocity_Inf_mean[iDim];
			}

			Pressure_Inf[iSpecies] = Pressure_Inf_mean * Gas_Composition[iSpecies];
			Density_Inf[iSpecies] = Pressure_Inf[iSpecies] / ((UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies]) * config->GetTemperature_FreeStream());
			Vel2 = 0.0;
			double Energy_el_Inf = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Inf[iSpecies][iDim] * Velocity_Inf[iSpecies][iDim];
			if (iSpecies < nDiatomics) {
				Energy_vib_Inf[iSpecies] = 0.0;
				Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*(GammaDiatomic-1.0)) + 1.0/2.0*Vel2 + Enthalpy_Formation[iSpecies] + Energy_vib_Inf[iSpecies] + Energy_el_Inf;
			}
			else {		
				Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*(GammaMonatomic-1.0)) + 1.0/2.0*Vel2 + Enthalpy_Formation[iSpecies] + Energy_el_Inf;
			}
		}

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			Density_Inlet[iSpecies] = Density_Inf[iSpecies];
			Density_Outlet[iSpecies] = Density_Inf[iSpecies];
			Energy_Inlet[iSpecies] = Energy_Inf[iSpecies];
			Energy_Outlet[iSpecies] = Energy_Inf[iSpecies];
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
			}
		}

		/*--- Restart the solution from file information ---*/
		if (!restart) {
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, Enthalpy_Formation, nDim, nVar, nSpecies,
						nMonatomics, nDiatomics, config);
		}
		else {
			unsigned long index;
			string text_line;

			/*--- Restart the solution from file information ---*/
			mesh_filename = config->GetSolution_FlowFileName();
			restart_file.open(mesh_filename.data(), ios::in);

			/*--- In case there is no file ---*/
			if (restart_file.fail()) {
				cout << "There is no flow restart file!!" << endl;
				cout << "Press any key to exit..." << endl;
				cin.get(); exit(1);
			}

			/*--- Read the restart file ---*/
			for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				getline(restart_file,text_line);
				istringstream point_line(text_line);
				point_line >> index;
				for (iVar = 0; iVar < nVar; iVar++)
					point_line >> Solution[iVar];
				node[iPoint] = new CPlasmaVariable(Solution, nDim, nVar, nSpecies, config);
			}

			/*--- Close the restart file ---*/
			restart_file.close();
		}
	}
	if (config->GetKind_GasModel() == AIR21) {

		double *Energy_Formation;
		nDim = geometry->GetnDim();
		nMonatomics = 19;
		nDiatomics = 0;
		nSpecies = nMonatomics + nDiatomics;
		nVar = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);
		node = new CVariable*[geometry->GetnPoint()];
		Atomic_Mass = new double[nSpecies];
		Gas_Constant = new double[nSpecies];

		double GammaDiatomic = config->GetGammaDiatomic();
		double GammaMonatomic = config->GetGammaMonatomic();
		double AtomicOxygen, AtomicNitrogen;

		AtomicOxygen = 2.6567625500E-26; AtomicNitrogen = 2.3258669700E-26;

		Atomic_Mass[0] = 	2*AtomicOxygen;
		Atomic_Mass[1] = 	3*AtomicOxygen;
		Atomic_Mass[2] = 	2*AtomicNitrogen;
		Atomic_Mass[3] = 	AtomicNitrogen+AtomicOxygen;
		Atomic_Mass[4] = 	2*AtomicNitrogen;
		Atomic_Mass[5] = 	2*AtomicNitrogen;
		Atomic_Mass[6] = 	2*AtomicNitrogen;
		Atomic_Mass[7] = 	2*AtomicNitrogen;
		Atomic_Mass[8] = 	2*AtomicNitrogen;
		Atomic_Mass[9] = 	2*AtomicOxygen-ELECTRON_MASS;
		Atomic_Mass[10] = 	2*AtomicNitrogen-ELECTRON_MASS;
		Atomic_Mass[11] = 	2*AtomicOxygen+ELECTRON_MASS;
		Atomic_Mass[12] = 	AtomicOxygen;
		Atomic_Mass[13] = 	AtomicNitrogen;
		Atomic_Mass[14] = 	AtomicOxygen;
		Atomic_Mass[15] = 	AtomicNitrogen;
		Atomic_Mass[16] = 	AtomicOxygen-ELECTRON_MASS;
		Atomic_Mass[17] = 	AtomicOxygen+ELECTRON_MASS;
		Atomic_Mass[18] = 	ELECTRON_MASS;

		Ru = UNIVERSAL_GAS_CONSTANT;
		double *Cv =  new double [nSpecies];


		/*--- Define  some auxiliar vector related with the residual ---*/
		Residual = new double[nVar];	Residual_Max = new double[nVar];
		Residual_i = new double[nVar];	Residual_j = new double[nVar];
		Res_Conv = new double[nVar];	Res_Visc = new double[nVar];	Res_Sour = new double[nVar];

		/*--- Define some auxiliar vector related with the solution ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];

		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector = new double[nDim];
		Vector_i = new double[nDim]; Vector_j = new double[nDim];

		/*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
		if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) {
			p1_Und_Lapl = new double [geometry->GetnPoint()]; p2_Und_Lapl = new double [geometry->GetnPoint()]; }

		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT) {
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar]; Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar]; Jacobian_j[iVar] = new double [nVar]; }
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Initialize_Jacobian_Structure(geometry, config);
			xsol = new double [geometry->GetnPoint()*nVar];
			rhs = new double [geometry->GetnPoint()*nVar];
		}

		/*--- Computation of gradients by least squares ---*/
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
				(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
			/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
			Smatrix = new double* [nDim];
			for (iDim = 0; iDim < nDim; iDim++)
				Smatrix[iDim] = new double [nDim];
			/*--- c vector := transpose(WA)*(Wb) ---*/
			cvector = new double* [nVar+nSpecies];
			for (iVar = 0; iVar < nVar+nSpecies; iVar++)
				cvector[iVar] = new double [nDim];
		}

		CSkinFriction = new double* [config->GetnMarker_All()];
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];

		CHeatTransfer = new double* [config->GetnMarker_All()];
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			CHeatTransfer[iMarker] = new double [geometry->nVertex[iMarker]];

		/*--- Flow infinity initialization stuff ---*/
		Density_Inf = new double [nSpecies];
		Pressure_Inf = new double [nSpecies];
		Velocity_Inf = new double*[nSpecies];
		Mach_Inf = 	new double [nSpecies];
		Energy_Inf = new double [nSpecies];
		Energy_vib_Inf = new double [nDiatomics];
		Temperature_Inf = new double [nSpecies];

		Energy_Formation = new double [nSpecies];
		Enthalpy_Formation = new double [nSpecies];

		Density_Inlet = new double [nSpecies];
		Velocity_Inlet  = new double*[nSpecies];
		Mach_Inlet 		= new double [nSpecies];
		Pressure_Inlet  = new double [nSpecies];
		Energy_Inlet = new double [nSpecies];
		Density_Outlet = new double [nSpecies];
		Pressure_Outlet = new double [nSpecies];
		Energy_Outlet = new double [nSpecies];
		Mach_Outlet = new double [nSpecies];
		Velocity_Outlet  = new double*[nSpecies];

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Velocity_Inf[iSpecies] = new double [nDim];
			Velocity_Inlet[iSpecies] = new double [nDim];
			Velocity_Outlet[iSpecies] = new double [nDim];
		}

#ifdef Air_AllSameValues

		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++ ){
			Density_Inf[iSpecies] = 1.2;
			Pressure_Inf[iSpecies] = 1E5;
			Velocity_Inf[iSpecies][0] = 100.0;
			for (iDim = 1; iDim < nDim; iDim ++)
				Velocity_Inf[iSpecies][iDim] = 0.0;
		}

		for ( iSpecies = 0; iSpecies < nDiatomics; iSpecies++ ){
			Gamma = GammaDiatomic;
			Gamma_Minus_One = Gamma -1.0;
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			Energy_vib_Inf[iSpecies] = 0.0;
			Energy_Formation[iSpecies] = 0.0;
			Enthalpy_Formation[iSpecies] = 0.0;
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];

			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*Gamma_Minus_One)+0.5*Vel2;
		}

		for ( iSpecies = nDiatomics; iSpecies < nSpecies; iSpecies++ ){
			Gamma = GammaMonatomic;
			Gamma_Minus_One = Gamma -1.0;
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			Energy_Formation[iSpecies] = 100.0;
			Enthalpy_Formation[iSpecies] = 100.0;
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];

			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*Gamma_Minus_One)+0.5*Vel2;
		}

		/*--- Flow at the inlet initialization stuff ---*/
		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Density_Inlet[iSpecies]  = Density_Inf[iSpecies];
			Mach_Inlet[iSpecies]     = Mach_Inf[iSpecies];
			Pressure_Inlet[iSpecies] = Pressure_Inf[iSpecies];
			Energy_Inlet[iSpecies]   = Energy_Inf[iSpecies];
			for (iDim = 0; iDim < nDim; iDim ++) {
				Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
			}
			Density_Outlet[iSpecies]  = Density_Inf[iSpecies];
			Mach_Outlet[iSpecies]     = Mach_Inf[iSpecies];
			Pressure_Outlet[iSpecies] = Pressure_Inf[iSpecies];
			Energy_Outlet[iSpecies]   = Energy_Inf[iSpecies];
		}
		cout << " nSpecies = " << nSpecies << endl;
		cout << " nDiatomics = " << nDiatomics << endl;
		cout << " Density_Inlet = " << Density_Inlet[0] << endl;
		cout << " Pressure_Inlet = " << Pressure_Inlet[0] << endl;
		cout << " Mach_Inlet = " << Mach_Inlet[0] << endl;

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, Enthalpy_Formation, nDim, nVar, nSpecies,
					nMonatomics, nDiatomics, config);

# endif

		//#ifdef Air_DifferentValues

		double total_density = 1.225;
		double total_pressure = 1E5;
		double fraction_excited = 1E-7;
		double fraction_pos_ion = 1E-7;
		double fraction_neg_ion = 1.5E-7;

		double *fraction = new double [nSpecies];
		fraction[0] = 0.21; fraction[1] = fraction_excited; fraction[2] = 0.79; fraction[3] = fraction_excited;
		fraction[4] = fraction_excited; fraction[5] = fraction_excited; fraction[6] = fraction_excited;
		fraction[7] = fraction_excited; fraction[8] = fraction_excited; fraction[9] = fraction_pos_ion;
		fraction[10] = fraction_pos_ion; fraction[11] = fraction_neg_ion;	fraction[12] = fraction_excited;
		fraction[13] = fraction_excited; fraction[14] = fraction_excited; fraction[15] = fraction_excited;
		fraction[16] = fraction_pos_ion; fraction[17] = fraction_neg_ion; fraction[18] = fraction_neg_ion;

		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++ ){
			Gas_Constant[iSpecies] = Ru/(Atomic_Mass[iSpecies] * AVOGAD_CONSTANT);
			Density_Inf[iSpecies] = fraction[iSpecies]*total_density*Atomic_Mass[iSpecies]/Atomic_Mass[0];
			Pressure_Inf[iSpecies] = fraction[iSpecies]*total_pressure;
			Temperature_Inf[iSpecies] = Pressure_Inf[iSpecies]/ (Gas_Constant[iSpecies] * Density_Inf[iSpecies]);
			Velocity_Inf[iSpecies][0] = 100.0;
			for (iDim = 1; iDim < nDim; iDim ++)
				Velocity_Inf[iSpecies][iDim] = 0.0;
		}

		for ( iSpecies = 0; iSpecies < nDiatomics; iSpecies++ ){
			Gamma = GammaDiatomic;
			Gamma_Minus_One = Gamma - 1.0;
			Cv[iSpecies]  = Gas_Constant[iSpecies]/(Gamma-1.0);
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			Energy_vib_Inf[iSpecies] = 0.0;
			Energy_Formation[iSpecies] = 0.0;
			Enthalpy_Formation[iSpecies] = 0.0;
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];

			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*Gamma_Minus_One)+0.5*Vel2;
		}

		for ( iSpecies = nDiatomics; iSpecies < nSpecies; iSpecies++ ){
			Gamma = GammaMonatomic;
			Gamma_Minus_One = Gamma -1.0;
			Cv[iSpecies]  = Gas_Constant[iSpecies]/(Gamma-1.0);
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(Gamma*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			Energy_Formation[iSpecies] = 100.0;
			Enthalpy_Formation[iSpecies] = 100.0;
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];

			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*Gamma_Minus_One)+0.5*Vel2;
		}

		/*--- Flow at the inlet initialization stuff ---*/
		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Density_Inlet[iSpecies]  = Density_Inf[iSpecies];
			Mach_Inlet[iSpecies]     = Mach_Inf[iSpecies];
			Pressure_Inlet[iSpecies] = Pressure_Inf[iSpecies];
			Energy_Inlet[iSpecies]   = Energy_Inf[iSpecies];
			for (iDim = 0; iDim < nDim; iDim ++) {
				Velocity_Inlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
				Velocity_Outlet[iSpecies][iDim] = Velocity_Inf[iSpecies][iDim];
			}
			Density_Outlet[iSpecies]  = Density_Inf[iSpecies];
			Mach_Outlet[iSpecies]     = Mach_Inf[iSpecies];
			Pressure_Outlet[iSpecies] = Pressure_Inf[iSpecies];
			Energy_Outlet[iSpecies]   = Energy_Inf[iSpecies];
		}
		cout << " nSpecies = " << nSpecies << endl;
		cout << " nDiatomics = " << nDiatomics << endl;
		cout << " Density_Inlet = " << Density_Inlet[0] << endl;
		cout << " Pressure_Inlet = " << Pressure_Inlet[0] << endl;
		cout << " Mach_Inlet = " << Mach_Inlet[0] << endl;

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
			cout << "iSpecies = " << iSpecies << " Density = " << Density_Inf[iSpecies] << " Mach = " << Mach_Inf[iSpecies] << " Temperature = " << Temperature_Inf[iSpecies] << " Pressure = " << Pressure_Inf[iSpecies] << endl;

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, Enthalpy_Formation, nDim, nVar, nSpecies,
					nMonatomics, nDiatomics, config);

		//# endif

	}
}

CPlasmaSolution::~CPlasmaSolution(void) {
	unsigned short iVar, iDim;

	delete [] Residual;		delete [] Residual_Max;
	delete [] Residual_i;	delete [] Residual_j;
	delete [] Residual_Chemistry;			delete [] Residual_ElecForce;
	delete [] Residual_MomentumExch;	delete[] Residual_EnergyExch;
	delete [] Res_Conv;		delete [] Res_Visc;		delete [] Res_Sour;
	delete [] Solution;		delete [] Solution_i; delete [] Solution_j;
	delete [] Vector;			delete [] Vector_i;		delete [] Vector_j;
	delete [] p1_Und_Lapl;  delete [] p2_Und_Lapl;
	delete [] xsol;			delete [] rhs;

	delete [] Density_Inf;
	delete [] Energy_Inf;
	delete [] Energy_vib_Inf;
	delete [] Pressure_Inf;
	delete [] Mach_Inf;
	delete [] Enthalpy_Formation;
	delete [] Molar_Mass;
	delete [] Gas_Composition;
	if (Mach2Vel_FreeStream != NULL) delete [] Mach2Vel_FreeStream;

	delete [] Energy_Inlet;		delete [] Energy_Outlet;
	delete [] Density_Inlet;	delete [] Density_Outlet;
	delete [] Pressure_Inlet;	delete [] Pressure_Outlet;
	delete [] Mach_Inlet;		delete [] Mach_Outlet;


	for (iVar = 0; iVar < nSpecies; iVar ++) {
		delete [] Velocity_Inf[iVar];
		delete [] Velocity_Inlet[iVar];
		delete [] Velocity_Outlet[iVar];
	}
	delete[] Velocity_Inf; 	delete [] Velocity_Inlet; delete [] Velocity_Outlet;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
		delete [] Jacobian_Chemistry[iVar];		
		delete [] Jacobian_ElecForce[iVar];
		delete [] Jacobian_MomentumExch[iVar];
		delete [] Jacobian_EnergyExch[iVar];		
	}
	delete [] Jacobian_i;							delete [] Jacobian_j;
	delete [] Jacobian_Chemistry;			delete [] Jacobian_ElecForce;
	delete [] Jacobian_MomentumExch;	delete [] Jacobian_EnergyExch;

	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;

	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;

}

void CPlasmaSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	double *Gas_Constant_MS = new double [nSpecies];
	bool viscous = (config->GetKind_ViscNumScheme_Plasma());

	if (config->GetKind_GasModel() == ARGON)
		for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
			Gas_Constant_MS[iSpecies] = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));



	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		/*--- Compute squared velocity, sound velocity, pressure, and enthalpy ---*/
		if (config->GetKind_GasModel() == ARGON) {
			node[iPoint]->SetVelocity2();
			node[iPoint]->SetSoundSpeed(Gamma);
			node[iPoint]->SetPressure(Gamma);
			node[iPoint]->SetTemperature(Gas_Constant_MS);
			node[iPoint]->SetEnthalpy();
			if(viscous) {
				node[iPoint]->SetLaminarViscosity(config);
				node[iPoint]->SetThermalCoeff(Gamma, Gas_Constant_MS);
			}	

		}

		if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == AIR21 ) {
			node[iPoint]->SetVelocity2();
			node[iPoint]->SetSoundSpeed(config->GetGammaMonatomic(), config->GetGammaDiatomic());
			node[iPoint]->SetPressure(config->GetGammaMonatomic(), config->GetGammaDiatomic(), geometry->node[iPoint]->GetCoord());
			node[iPoint]->SetTemperature_TR(Molar_Mass, config->GetGammaMonatomic(), config->GetGammaDiatomic());
			node[iPoint]->SetEnthalpy();
		}


		/*--- Initialize the convective and viscous residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResSour_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit)
			node[iPoint]->Set_ResVisc_Zero();
	}


	if (config->GetKind_SourNumScheme_Plasma() !=NONE)
		SetSolution_Gradient_GG(geometry);

	/*--- Initialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT)
		Jacobian.SetValZero();

}

void CPlasmaSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Area, dV, Mean_SoundSpeed, Mean_ProjVel, Lambda, Lambda_1, Lambda_2, Local_Delta_Time, Global_Delta_Time = 1E6;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker, iSpecies, loc;
	double Lambda_iSpecies, Mean_LaminarVisc, Mean_Density;

	double *Gas_Constant_MS = new double [nSpecies];

	double Lambda_Visc_iSpecies,K_v, Local_Delta_Time_Visc;

	K_v = 0.25;
	bool viscous = (config->GetKind_ViscNumScheme_Plasma());
	bool MultipleTimeSteps = (config->MultipleTimeSteps());

	if (config->GetKind_GasModel() == ARGON) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
			Gas_Constant_MS[iSpecies] = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));

		/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				node[iPoint]->SetMax_Lambda_Inv(0.0, iSpecies);
				node[iPoint]->SetMax_Lambda_Visc(0.0, iSpecies);
			}
			node[iPoint]->SetVelocity2();
			node[iPoint]->SetPressure(Gamma);
			node[iPoint]->SetSoundSpeed(Gamma);
			if (viscous) {
				node[iPoint]->SetTemperature(Gas_Constant_MS);
				node[iPoint]->SetLaminarViscosity(config);
			}
		}
		/*--- Loop interior edges ---*/
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			Normal = geometry->edge[iEdge]->GetNormal();
			Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

			/*--- Mean Values ---*/
			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
				Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal,iSpecies) + node[jPoint]->GetProjVel(Normal,iSpecies));
				Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed(iSpecies) + node[jPoint]->GetSoundSpeed(iSpecies));
				Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
				/*--- Inviscid contribution ---*/
				node[iPoint]->AddMax_Lambda_Inv(Lambda,iSpecies);
				node[jPoint]->AddMax_Lambda_Inv(Lambda,iSpecies);

				if (viscous) {
					loc = (nDim+2)*iSpecies;
					Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity(iSpecies) + node[jPoint]->GetLaminarViscosity(iSpecies));
					Mean_Density     = 0.5*(node[iPoint]->GetSolution(loc + 0) + node[jPoint]->GetSolution(loc + 0));

					Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
					Lambda_2 = Gamma*Mean_LaminarVisc/Prandtl_Lam;
					Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

					node[iPoint]->AddMax_Lambda_Visc(Lambda,iSpecies);
					node[jPoint]->AddMax_Lambda_Visc(Lambda,iSpecies);

				}
			}
		}

		/*--- Loop boundary edges ---*/
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

				/*--- Point identification, Normal vector and area ---*/
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

				for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
					/*--- Mean Values ---*/
					Mean_ProjVel = node[iPoint]->GetProjVel(Normal,iSpecies);
					Mean_SoundSpeed = node[iPoint]->GetSoundSpeed(iSpecies);
					/*--- Inviscid contribution ---*/
					Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
					node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
					if (viscous) {
						loc = (nDim+2)*iSpecies;
						Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity(iSpecies);
						Mean_Density     = node[iPoint]->GetSolution(loc + 0);
						Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
						Lambda_2 = Gamma*Mean_LaminarVisc/Prandtl_Lam;
						Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
						node[iPoint]->AddMax_Lambda_Visc(Lambda,iSpecies);
					}
				}
			}
		}
		/* 	double u, c, dx; u = 4800; c = 732; dx = 0.04/810; Local_Delta_Time = config->GetCFL(iMesh)*dx/(u+c); */
		if (!MultipleTimeSteps) {
			/*--- Each element uses their own speed ---*/
			for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
				dV = geometry->node[iPoint]->GetVolume();
				Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(0);
				for (iSpecies = 1; iSpecies < nSpecies; iSpecies++)
					Lambda_iSpecies = max(Lambda_iSpecies, node[iPoint]->GetMax_Lambda_Inv(iSpecies));

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
			/*--- Each element uses their own speed ---*/
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

	if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == AIR21) {
		/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++ )
				node[iPoint]->SetMax_Lambda_Inv(0.0, iSpecies);

			node[iPoint]->SetVelocity2();
			node[iPoint]->SetSoundSpeed(config->GetGammaMonatomic(), config->GetGammaDiatomic());
		}

		/*--- Loop interior edges ---*/
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			Normal = geometry->edge[iEdge]->GetNormal();
			Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

			/*--- Mean Values ---*/
			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
				Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal,iSpecies, nDiatomics) + node[jPoint]->GetProjVel(Normal,iSpecies,nDiatomics));
				Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed(iSpecies) + node[jPoint]->GetSoundSpeed(iSpecies));
				Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;

				/*--- Inviscid contribution ---*/
				node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
				node[jPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
			}
		}

		/*--- Loop boundary edges ---*/
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

				/*--- Point identification, Normal vector and area ---*/
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

				for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {

					/*--- Mean Values ---*/
					Mean_ProjVel = node[iPoint]->GetProjVel(Normal,iSpecies,nDiatomics);
					Mean_SoundSpeed = node[iPoint]->GetSoundSpeed(iSpecies);

					/*--- Inviscid contribution ---*/
					Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
					node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
				}
			}
		}

		/*--- Each element uses their own speed ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			dV = geometry->node[iPoint]->GetVolume();

			Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(0);
			for (iSpecies = 1; iSpecies < nSpecies; iSpecies++)
				Lambda_iSpecies = max(Lambda_iSpecies, node[iPoint]->GetMax_Lambda_Inv(iSpecies));

			Local_Delta_Time = config->GetCFL(iMesh)*dV / Lambda_iSpecies;
			if (config->GetUnsteady_Simulation() == TIME_STEPPING) Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
			else node[iPoint]->SetDelta_Time(Local_Delta_Time);

		}

		/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
		if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
			for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
				node[iPoint]->SetDelta_Time(Global_Delta_Time);
		}
	}
}

void CPlasmaSolution::Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CPlasmaSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Upwind() == ROE_2ND) && (iMesh == MESH_0));

	//SetSolution_Gradient_GG(geometry);

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
				(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) SetSolution_Gradient_LS(geometry, config);
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
			Gradient_i = node[iPoint]->GetGradient();
			Gradient_j = node[jPoint]->GetGradient();
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_j = node[jPoint]->GetLimiter();
				Limiter_i = node[iPoint]->GetLimiter();
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
			/*--- Set conservative variables with reconstruction ---*/
			solver->SetConservative(Solution_i, Solution_j);
		}

		/*--- Compute the residual ---*/

		solver->SetResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

		/*		if (iPoint < 5) {
			cout << " iPoint = " << iPoint << endl;
			for (iVar = 0; iVar < nVar; iVar ++)
				cout << "  " << iVar << " " << Res_Conv[iVar] << endl;
		}
		 */

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
	double **Gradient;
	double *Coord_0 = NULL, *Coord_1= NULL;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);


	if (config->GetKind_GasModel() == ARGON || config->GetKind_GasModel() == AIR21) {

		for (iVar = 0; iVar < nVar; iVar++)
			Residual[iVar] = 0;

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

			if (config->GetElectricSolver()) {
				Gradient = solution_container[ELEC_SOL]->node[iPoint]->GetGradient();
				/*--- Set gradient of phi from electrostatic potential equation ---*/
				solver->SetConsVarGradient(Gradient);
			}
			Coord_0 = geometry->node[iPoint]->GetCoord();
			Coord_1 = geometry->node[iPoint]->GetCoord();

			solver->SetCoord(Coord_0, Coord_1);

			/*--- Set solution  ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

			/*--- Set control volume ---*/
			solver->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Compute Residual ---*/
			solver->SetResidual(Residual,Jacobian_i, config);
			solution_container[PLASMA_SOL]->node[iPoint]->SetSource(Residual);


			/*--- Add Residual ---*/
			node[iPoint]->SubtractRes_Sour(Residual);

			if (implicit) {
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				//		solution_container[PLASMA_SOL]->node[iPoint]->SetSourceJacobian(Jacobian_i);

			}
		}
	}

	if (config->GetKind_GasModel() == AIR7) {
		double *Temperature_tr_i;

		Temperature_tr_i = new double [nSpecies];

		for (iVar = 0; iVar < nVar; iVar++) {
			Residual[iVar] = 0.0;
			Residual_Chemistry[iVar] = 0.0;
			Residual_MomentumExch[iVar] = 0.0;
			Residual_ElecForce[iVar] = 0.0;
			Residual_EnergyExch[iVar] = 0.0;
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_Chemistry[iVar][jVar] = 0.0;
				Jacobian_MomentumExch[iVar][jVar] = 0.0;
			}
		}

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
				Temperature_tr_i[iSpecies] = node[iPoint]->GetTemperature_TR(iSpecies);
			}

			Gradient = solution_container[ELEC_SOL]->node[iPoint]->GetGradient();

			/*--- Set gradient of phi from electrostatic potential equation ---*/
			solver->SetConsVarGradient(Gradient);

			Coord_0 = geometry->node[iPoint]->GetCoord();
			Coord_1 = geometry->node[iPoint]->GetCoord();

			solver->SetCoord(Coord_0, Coord_1);

			/*--- Set solution  ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
			Gradient = solution_container[ELEC_SOL]->node[iPoint]->GetGradient();
			solver->SetConsVarGradient(Gradient);

			/*--- Set control volume ---*/
			solver->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Set temperature ---*/
			solver->SetTemperature_TR(Temperature_tr_i, Temperature_tr_i);

			/*--- Compute Residual ---*/
			solver->SetResidual_Chemistry(Residual_Chemistry, config);
			solver->SetJacobian_Chemistry(Jacobian_Chemistry, config);
			node[iPoint]->SubtractRes_Sour(Residual_Chemistry);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_Chemistry);

			/*
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << Jacobian_ElecForce[iVar][jVar];
			}
			cout << endl;
		}
		cout <<endl;

		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << Jacobian_Chemistry[iVar][jVar];
			}
			cout << endl;
		}
		cin.get(); */

			/*		solver->SetResidual_ElecForce(Residual_ElecForce, config);
		 solver->SetJacobian_ElecForce(Jacobian_ElecForce, config);
		 node[iPoint]->SubtractRes_Sour(Residual_ElecForce);
		 if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ElecForce);*/

			/*		solver->SetResidual_MomentumExch(Residual_MomentumExch, config);
		solver->SetJacobian_MomentumExch(Jacobian_MomentumExch, config);
		node[iPoint]->SubtractRes_Sour(Residual_MomentumExch);
		if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_MomentumExch);*/

			/*		solver->SetResidual_EnergyExch(Residual_EnergyExch, Residual_ElecForce, config);
		solver->SetJacobian_EnergyExch(Jacobian_EnergyExch, config);
		node[iPoint]->SubtractRes_Sour(Residual_EnergyExch);
		if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_EnergyExch); */

		}
		delete[] Temperature_tr_i;
	}
}

void CPlasmaSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge;
	double Pressure_i = 0, Pressure_j = 0;
	unsigned short iVar;
	bool boundary_i, boundary_j;
	double *Diff = new double[nVar];

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetUnd_LaplZero();

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		Pressure_i = node[iPoint]->GetPressure();
		Pressure_j = node[jPoint]->GetPressure();

		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
		Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1)+Pressure_i) - (node[jPoint]->GetSolution(nVar-1)+Pressure_j);

		boundary_i = geometry->node[iPoint]->GetBoundary_Physical();
		boundary_j = geometry->node[jPoint]->GetBoundary_Physical();

		/*--- Both points inside Omega ---*/
		if (!boundary_i && !boundary_j) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}

		/*--- iPoint inside Omega, jPoint on the boundary ---*/
		if (!boundary_i && boundary_j)
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);

		/*--- jPoint inside Omega, iPoint on the boundary ---*/
		if (boundary_i && !boundary_j)
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

		/*--- Both points on the boundary ---*/
		if (boundary_i && boundary_j) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}
	}

	delete [] Diff;
}

void CPlasmaSolution::SetSpectral_Radius(CGeometry *geometry, CConfig *config) {
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;
	double *Normal, Area, ProjVel, Lambda, SoundSpeed_i, SoundSpeed_j, SoundSpeed;

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetLambda(0.0);

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		SoundSpeed_i = node[iPoint]->GetSoundSpeed();
		SoundSpeed_j = node[jPoint]->GetSoundSpeed();

		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

		/*--- Inviscid contribution to the Point i ---*/
		ProjVel = node[iPoint]->GetProjVel(Normal);
		Lambda = 0.5*(fabs(ProjVel) + SoundSpeed_i*Area);
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);

		/*--- Inviscid contribution to the Point j ---*/
		ProjVel = node[jPoint]->GetProjVel(Normal);
		Lambda = 0.5*(fabs(ProjVel) + SoundSpeed_j*Area);
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
	}

	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {
				SoundSpeed = node[iPoint]->GetSoundSpeed();
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
				/*--- Inviscid contribution to the iPoint ---*/
				ProjVel = node[iPoint]->GetProjVel(Normal);
				Lambda = (fabs(ProjVel) + SoundSpeed*Area);
				node[iPoint]->AddLambda(Lambda);
			}
		}
}

void CPlasmaSolution::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {
	unsigned long iEdge, iPoint, jPoint;
	double Pressure_i, Pressure_j;

	/*--- Reset variables ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		p1_Und_Lapl[iPoint] = 0.0;
		p2_Und_Lapl[iPoint] = 0.0;
	}

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		Pressure_i = node[iPoint]->GetPressure();
		Pressure_j = node[jPoint]->GetPressure();

		if (geometry->node[iPoint]->GetDomain()) p1_Und_Lapl[iPoint] += (Pressure_j - Pressure_i);
		if (geometry->node[jPoint]->GetDomain()) p1_Und_Lapl[jPoint] += (Pressure_i - Pressure_j);

		if (geometry->node[iPoint]->GetDomain()) p2_Und_Lapl[iPoint] += (Pressure_i + Pressure_j);
		if (geometry->node[jPoint]->GetDomain()) p2_Und_Lapl[jPoint] += (Pressure_i + Pressure_j);

	}

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetSensor(fabs(p1_Und_Lapl[iPoint])/p2_Und_Lapl[iPoint]);
}

void CPlasmaSolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container,
		CConfig *config, unsigned short iRKStep) {
	double *Residual, *TruncationError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;
		TruncationError = node[iPoint]->GetTruncationError();
		Residual = node[iPoint]->GetResidual();
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -(Residual[iVar]+TruncationError[iVar])*Delta*RK_AlphaCoeff);
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
	double *Residual, *TruncationError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;

	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;
		TruncationError = node[iPoint]->GetTruncationError();
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
	double Delta, Res, *local_ResConv, *local_ResVisc, *local_ResSour, *local_TruncationError, Vol;
	double *Species_Delta;
	bool MultipleTimeSteps = (config->MultipleTimeSteps());

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);

	Species_Delta = new double [nSpecies];

	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		local_TruncationError = node[iPoint]->GetTruncationError();
		local_ResConv = node[iPoint]->GetResConv();
		local_ResVisc = node[iPoint]->GetResVisc();
		local_ResSour = node[iPoint]->GetResSour();
		Vol = geometry->node[iPoint]->GetVolume();

		if (!MultipleTimeSteps) {
			/*--- Modify matrix diagonal to assure diagonal dominance ---*/
			Delta = geometry->node[iPoint]->GetVolume()/node[iPoint]->GetDelta_Time();
			Jacobian.AddVal2Diag(iPoint,Delta);
		} else {
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
				Species_Delta[iSpecies] = geometry->node[iPoint]->GetVolume()/node[iPoint]->GetDelta_Time(iSpecies);
			}
			Jacobian.AddVal2Diag(iPoint, Species_Delta, nDim);
		}

		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
			Res = local_ResConv[iVar]+local_ResVisc[iVar]+local_ResSour[iVar];
			rhs[total_index] = -(Res+local_TruncationError[iVar]);
			xsol[total_index] = 0.0;
			AddRes_Max( iVar, Res*Res*Vol );
		}

	}

	/*--- Solve the system ---*/
	if(config->GetKind_GasModel() == AIR7) {
		if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) Jacobian.SGSSolution(rhs, xsol, config->GetLinear_Solver_Error(),
				config->GetLinear_Solver_Iter(), true, geometry, config);
		if (config->GetKind_Linear_Solver() == BCGSTAB) Jacobian.BCGSTABSolution(rhs, xsol, config->GetLinear_Solver_Error(),
				config->GetLinear_Solver_Iter(), true, geometry, config);
		if (config->GetKind_Linear_Solver() == LU_SGS) Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
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
		if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) Jacobian.SGSSolution(rhs, xsol, 1e-9, 100, true, geometry, config);
		if (config->GetKind_Linear_Solver() == LU_SGS) Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			for (iVar = 0; iVar < nVar; iVar++) {
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
			loc = (nDim+3)*iSpecies;
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
					loc = (nDim+3)*iSpecies;
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
				PrimVar_j[loc + 0] = node[jPoint]->GetTemperature(iSpecies);
				for (iDim = 0; iDim < nDim; iDim++)
					PrimVar_j[loc + iDim+1] = node[jPoint]->GetVelocity(iDim,iSpecies);
				PrimVar_j[loc + nDim+1] = node[jPoint]->GetPressure(iSpecies);
				PrimVar_j[loc + nDim+2] = node[jPoint]->GetDensity(iSpecies);
			}
			weight = 0.0;
			//				if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
			for (iDim = 0; iDim < nDim; iDim++)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
			//				}
			//				else weight = 1.0;

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
				for (jDim = 0; jDim < nDim; jDim++)
					product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
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
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
					if ( iSpecies < nDiatomics ) a2 = config->GetGammaDiatomic()-1.0;
					else  a2 = config->GetGammaMonatomic()-1.0;
					phi = 0.5*a2*node[iPoint]->GetVelocity2(iSpecies);

					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[loc + 0][loc + iVar] = 0.0;
						Jacobian_i[loc + nDim+1][loc + iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
					}
				}

				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
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
						loc = (nDim+2)*iSpecies;
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


void CPlasmaSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
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
					loc = (nDim+2)*iSpecies;
					node[iPoint]->SetVelocity_Old(Vector, iSpecies);
					node[iPoint]->SetVel_ResConv_Zero(iSpecies);
					node[iPoint]->SetVel_ResVisc_Zero(iSpecies);
					node[iPoint]->SetVelTruncationErrorZero(iSpecies);

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

				//CHECK THIS ONE
				//				Point_Normal = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();

				/*--- Vector --> Velocity_corrected ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;


				for (iVar = 0; iVar < nVar; iVar ++)
					Res_Visc[iVar] = 0.0;

				/*--- Set the residual, truncation error and velocity value ---*/
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					node[iPoint]->SetVelocity_Old(Vector);
					node[iPoint]->SetVel_ResConv_Zero(iSpecies);
					node[iPoint]->SetVel_ResVisc_Zero(iSpecies);
					node[iPoint]->SetVelTruncationErrorZero(iSpecies);
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
					loc = (nDim+2)*iSpecies;
					Temperature = node[Point_Normal]->GetTemperature(iSpecies);
					Temperature_Gradient = (Twall - Temperature)/dist_ij;
					Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
					gas_Constant = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));
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
						loc = (nDim+2)* iSpecies;
						sq_vel = 0.0;
						for (iDim = 0; iDim< nDim; iDim ++)
							sq_vel += (U_i[loc + iDim+1]/U_i[loc + 0])*(U_i[loc + iDim+1]/U_i[loc + 0]);

						Density = node[iPoint]->GetDensity(iSpecies);
						Pressure = node[iPoint]->GetPressure(iSpecies);
						Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
						gas_Constant = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));
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
					node[iPoint]->SetVelocity_Old(Vector);
					node[iPoint]->SetVel_ResConv_Zero(iSpecies);
					node[iPoint]->SetVel_ResVisc_Zero(iSpecies);
					node[iPoint]->SetVelTruncationErrorZero(iSpecies);
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
					loc = (nDim+2)*iSpecies;
					Temperature = node[Point_Normal]->GetTemperature(iSpecies);
					Temperature_Gradient = (Twall - Temperature)/dist_ij;
					Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
					gas_Constant = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));
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
						loc = (nDim+2)* iSpecies;
						sq_vel = 0.0;
						for (iDim = 0; iDim< nDim; iDim ++)
							sq_vel += (U_i[loc + iDim+1]/U_i[loc + 0])*(U_i[loc + iDim+1]/U_i[loc + 0]);

						Density = node[iPoint]->GetDensity(iSpecies);
						Pressure = node[iPoint]->GetPressure(iSpecies);
						Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
						gas_Constant = UNIVERSAL_GAS_CONSTANT/(AVOGAD_CONSTANT*config->GetParticle_Mass(iSpecies));
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
					/*
					iSpecies = nSpecies - 1;
					loc = (nDim+2)* iSpecies;
					for (iVar = loc + 1; iVar <= loc + nDim; iVar++) {
						total_index = iPoint*nVar+iVar;
						Jacobian.DeleteValsRowi(total_index);
					}*/
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
				loc = iSpecies * (nDim+2);
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

				/*		rho[iSpecies] = U_inlet[loc + 0];
				rhoE[iSpecies] = U_inlet[loc + nDim + 1];
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iSpecies][iDim] = U_inlet[loc + iDim+1]/U_inlet[loc + 0];
					sq_vel += velocity[iSpecies][iDim]*velocity[iSpecies][iDim];
					vn[iSpecies] += velocity[iSpecies][iDim]*kappa[iDim]*dS;
				}
				 */
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

				/*	rho[iSpecies]  = U_outlet[loc+0];
				rhoE[iSpecies] = U_outlet[loc + nDim + 1];
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iSpecies][iDim] = U_outlet[loc + iDim+1]/U_outlet[loc + 0];
					sq_vel += velocity[iSpecies][iDim]*velocity[iSpecies][iDim];
					vn[iSpecies] += velocity[iSpecies][iDim]*UnitaryNormal[iDim]*Area;
				}
				 */
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
				c[iSpecies] = sqrt(Gamma*Gamma_Minus_One*(rhoE[iSpecies]/rho[iSpecies] - 0.5*sq_vel));
				pressure[iSpecies] = (c[iSpecies] * c[iSpecies] * rho[iSpecies]) / Gamma;
				enthalpy[iSpecies] = (rhoE[iSpecies] + pressure[iSpecies]) / rho[iSpecies];
			}

			solver->GetInviscidProjFlux(rho, velocity, pressure, enthalpy, UnitaryNormal, Residual);

			for	(iVar = 0; iVar < nVar; iVar++)
				Residual[iVar] = Residual[iVar]*Area;

			node[iPoint]->AddRes_Conv(Residual);

#ifdef Argon1D
			node[iPoint]->Set_ResConv_Zero();
			node[iPoint]->Set_ResVisc_Zero();
			node[iPoint]->Set_ResSour_Zero();
#endif
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

	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double 	*UnitaryNormal, *U_domain;
	double *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	UnitaryNormal = new double[nDim];
	U_domain = new double[nVar];

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);

			/*--- Compute the projected residual ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_domain);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);

			node[iPoint]->Set_ResConv_Zero();
			node[iPoint]->Set_ResSour_Zero();
			node[iPoint]->Set_ResVisc_Zero();
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

		}
	}
	delete [] UnitaryNormal;
	delete [] U_domain;
}

#ifdef Oldfar
void CPlasmaSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim, iSpecies;
	unsigned short loc = 0;
	double *U_wall, *U_infty;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	U_wall = new double[nVar]; U_infty = new double[nVar];

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_wall[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the infinity ---*/
			for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				if (nDim == 2) {
					U_infty[loc + 0] = GetDensity_Inf(iSpecies);
					U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iSpecies);
					U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iSpecies);
					U_infty[loc + 3] = GetDensity_Energy_Inf(iSpecies);
					if ( iSpecies < nDiatomics ) U_infty[loc + 4] = GetDensity_Energy_vib_Inf(iSpecies);
				}
				if (nDim == 3) {
					U_infty[loc + 0] = GetDensity_Inf(iSpecies);
					U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iSpecies);
					U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iSpecies);
					U_infty[loc + 3] = GetDensity_Velocity_Inf(2,iSpecies);
					U_infty[loc + 4] = GetDensity_Energy_Inf(iSpecies);
					if ( iSpecies < nDiatomics ) U_infty[loc + 5] = GetDensity_Energy_vib_Inf(iSpecies);
				}
			}

			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];

			solver->SetNormal(Vector);
			solver->SetConservative(U_wall, U_infty);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

		}
	}
}
#endif

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
				U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iSpecies);
				U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iSpecies);
				U_infty[loc + 3] = GetDensity_Energy_Inf(iSpecies);
				if (nDim == 3) {
					U_infty[loc + 3] = GetDensity_Velocity_Inf(2,iSpecies);
					U_infty[loc + 4] = GetDensity_Energy_Inf(iSpecies);
				}
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
}

/*!
 * \method BC_Sym_Plane
 * \brief Symmetry Boundary Condition
 * \author A. Lonkar
 */

void CPlasmaSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar, iSpecies, loc;
	double Pressure, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Compute the projected residual ---*/
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
					if (config->GetKind_GasModel() == ARGON) {
						a2 = Gamma - 1.0;
					}
					if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == AIR21 ) {
						if ( iSpecies < nDiatomics ) a2 = config->GetGammaDiatomic()-1.0;
						else  a2 = config->GetGammaMonatomic()-1.0;
					}
					phi = 0.5*a2*node[iPoint]->GetVelocity2(iSpecies);

					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

					for (iVar =  0; iVar <  nVar; iVar++) {
						Jacobian_i[loc + 0][iVar] = 0.0;
						Jacobian_i[loc + nDim+1][iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
					}
				}

				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}

		}
	}
}

void CPlasmaSolution::MPI_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config,
		unsigned short val_mesh) {

#ifndef NO_MPI
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint;
	double *Conserv_Var, **Conserv_Grad = NULL;

	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

			short SendRecv = config->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry->nVertex[iMarker];

			/*--- Send information  ---*/
			if (SendRecv > 0) {

				/*--- Upwind scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {

					unsigned long nBuffer = geometry->nVertex[iMarker]*nVar;
					int send_to = SendRecv-1;

					double *Buffer_Send_U = new double[nBuffer];
					double *Buffer_Send_Ux = new double[nBuffer];
					double *Buffer_Send_Uy = new double[nBuffer];
					double *Buffer_Send_Uz = new double[nBuffer];
					//		double *Buffer_Send_Limit = new double[nBuffer];

					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						Conserv_Var = node[iPoint]->GetSolution();
						if (val_mesh == MESH_0) Conserv_Grad = node[iPoint]->GetGradient();
						/*				if ((val_mesh == MESH_0) && (config->GetKind_SlopeLimit() != NONE)) // make this plasma
						 Grad_Limit = node[iPoint]->GetLimiter();*/

						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_U[iVar*nVertex+iVertex] = Conserv_Var[iVar];
							if (val_mesh == MESH_0) {
								Buffer_Send_Ux[iVar*nVertex+iVertex] = Conserv_Grad[iVar][0];
								Buffer_Send_Uy[iVar*nVertex+iVertex] = Conserv_Grad[iVar][1];
								if (nDim == 3) Buffer_Send_Uz[iVar*nVertex+iVertex] = Conserv_Grad[iVar][2];
								/*if (config->GetKind_SlopeLimit() != NONE)
								 Buffer_Send_Limit[iVar*nVertex+iVertex] = Grad_Limit[iVar];*/
							}
						}
					}

					MPI::COMM_WORLD.Bsend(Buffer_Send_U,nBuffer,MPI::DOUBLE,send_to, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Bsend(Buffer_Send_Ux,nBuffer,MPI::DOUBLE,send_to, 1);
						MPI::COMM_WORLD.Bsend(Buffer_Send_Uy,nBuffer,MPI::DOUBLE,send_to, 2);
						if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_Uz,nBuffer,MPI::DOUBLE,send_to, 3);
						//				if (config->GetKind_SlopeLimit() != NONE) MPI::COMM_WORLD.Bsend(Buffer_Send_Limit,nBuffer,MPI::DOUBLE,send_to, 4);
					}

					delete [] Buffer_Send_U;
					delete [] Buffer_Send_Ux;
					delete [] Buffer_Send_Uy;
					delete [] Buffer_Send_Uz;
					//		delete [] Buffer_Send_Limit;

				}
			}

			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				/*--- Upwind scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {

					unsigned long nBuffer = geometry->nVertex[iMarker]*nVar;
					int receive_from = abs(SendRecv)-1;

					double *Buffer_Receive_U = new double[nBuffer];
					double *Buffer_Receive_Ux = new double[nBuffer];
					double *Buffer_Receive_Uy = new double[nBuffer];
					double *Buffer_Receive_Uz = new double[nBuffer];
					//		double *Buffer_Receive_Limit = new double[nBuffer];

					MPI::COMM_WORLD.Recv(Buffer_Receive_U,nBuffer,MPI::DOUBLE,receive_from, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Recv(Buffer_Receive_Ux,nBuffer,MPI::DOUBLE,receive_from, 1);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Uy,nBuffer,MPI::DOUBLE,receive_from, 2);
						if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_Uz,nBuffer,MPI::DOUBLE,receive_from, 3);
						//				if (config->GetKind_SlopeLimit() != NONE) MPI::COMM_WORLD.Recv(Buffer_Receive_Limit,nBuffer,MPI::DOUBLE,receive_from, 4); // change to plasma
					}

					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
							if (val_mesh == MESH_0) {
								node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Ux[iVar*nVertex+iVertex]);
								node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Uy[iVar*nVertex+iVertex]);
								if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Uz[iVar*nVertex+iVertex]);
								//	if (config->GetKind_SlopeLimit() != NONE) node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVar*nVertex+iVertex]);  // change to plasma
							}

						}
					}

					delete [] Buffer_Receive_U;
					delete [] Buffer_Receive_Ux;
					delete [] Buffer_Receive_Uy;
					delete [] Buffer_Receive_Uz;
					//			delete [] Buffer_Receive_Limit;
				}
			}
		}
#endif
}
