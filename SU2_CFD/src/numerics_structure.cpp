/*!
 * \file numerics_structure.cpp
 * \brief This file contains all the numerical methods.
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

#include "../include/numerics_structure.hpp"

CNumerics::CNumerics(void) { }

CNumerics::CNumerics(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) {
	nDim = val_nDim;
	nVar = val_nVar;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	Gas_Constant = config->GetGas_Constant();

	UnitaryNormal = new double [nDim];
	Normal = new double [nDim];
	Flux_Tensor = new double* [nVar];
	for (unsigned short iVar = 0; iVar < (nVar); iVar++)
		Flux_Tensor[iVar] = new double [nDim];

	tau = new double* [nDim];
	delta = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		tau[iDim] = new double [nDim];
		delta[iDim] = new double [nDim];
	}

	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		for (unsigned short jDim = 0; jDim < nDim; jDim++) {
			if (iDim==jDim) delta[iDim][jDim]=1.0;
			else delta[iDim][jDim]=0.0;
		}
	}

	U_n = new double [nVar];
	U_nM1 = new double [nVar];
	U_nP1 = new double [nVar];

	Proj_Flux_Tensor = new double [nVar];

}
/* Class overloaded to include multiple fluid equations for plasma */
CNumerics::CNumerics(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nFluids, CConfig *config) {
	nDim		= val_nDim;
	nVar		= val_nVar;
	nSpecies	= val_nSpecies;
	nFluids		= val_nFluids;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;


	UnitaryNormal = new double [nDim];
	Normal  = new double [nDim];
	Flux_Tensor = new double* [nVar];
	for (unsigned short iVar = 0; iVar < (nVar); iVar++)
		Flux_Tensor[iVar] = new double [nDim];

	tau = new double* [nDim];
	delta = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		tau[iDim] = new double [nDim];
		delta[iDim] = new double [nDim];
	}

	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		for (unsigned short jDim = 0; jDim < nDim; jDim++) {
			if (iDim==jDim) delta[iDim][jDim]=1.0;
			else delta[iDim][jDim]=0.0;
		}
	}

	U_n = new double [nVar];
	U_nM1 = new double [nVar];
	U_nP1 = new double [nVar];

	Proj_Flux_Tensor = new double [nVar];

}

/* Class overloaded to include multiple fluid equations for plasma */
CNumerics::CNumerics(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) {
	nDim		= val_nDim;
	nVar		= val_nVar;
	nSpecies	= val_nSpecies;
	nMonatomics = val_nMonatomics;
	nDiatomics  = val_nDiatomics;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	UnitaryNormal = new double [nDim];
	Normal  = new double [nDim];
	cout << " in constructor nDim = " << nDim << endl;
	Flux_Tensor = new double* [nVar];
	for (unsigned short iVar = 0; iVar < (nVar); iVar++)
		Flux_Tensor[iVar] = new double [nDim];

	tau = new double* [nDim];
	delta = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		tau[iDim] = new double [nDim];
		delta[iDim] = new double [nDim];
	}

	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		for (unsigned short jDim = 0; jDim < nDim; jDim++) {
			if (iDim==jDim) delta[iDim][jDim]=1.0;
			else delta[iDim][jDim]=0.0;
		}
	}

	U_n = new double [nVar];
	U_nM1 = new double [nVar];
	U_nP1 = new double [nVar];

	Proj_Flux_Tensor = new double [nVar];

}

CNumerics::~CNumerics(void) {

	delete [] UnitaryNormal;

	delete [] U_n;
	delete [] U_nM1;
	delete [] U_nP1;

	// visc	
	delete [] Proj_Flux_Tensor;

	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		delete [] Flux_Tensor[iVar];
	}
	delete [] Flux_Tensor;

	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		delete [] tau[iDim];
		delete [] delta[iDim];
	}
	delete [] tau;
	delete [] delta;	
	delete [] Normal;
	//cout << "~CNumerics()"<<endl;
}

void CNumerics::GetInviscidFlux(double val_density, double *val_velocity, 
		double val_pressure, double val_enthalpy) {
	if (nDim == 3) {
		Flux_Tensor[0][0] = val_density*val_velocity[0]; 
		Flux_Tensor[1][0] = Flux_Tensor[0][0]*val_velocity[0]+val_pressure; 
		Flux_Tensor[2][0] = Flux_Tensor[0][0]*val_velocity[1];
		Flux_Tensor[3][0] = Flux_Tensor[0][0]*val_velocity[2];
		Flux_Tensor[4][0] = Flux_Tensor[0][0]*val_enthalpy;

		Flux_Tensor[0][1] = val_density*val_velocity[1];
		Flux_Tensor[1][1] = Flux_Tensor[0][1]*val_velocity[0];
		Flux_Tensor[2][1] = Flux_Tensor[0][1]*val_velocity[1]+val_pressure; 
		Flux_Tensor[3][1] = Flux_Tensor[0][1]*val_velocity[2];
		Flux_Tensor[4][1] = Flux_Tensor[0][1]*val_enthalpy;

		Flux_Tensor[0][2] = val_density*val_velocity[2];
		Flux_Tensor[1][2] = Flux_Tensor[0][2]*val_velocity[0];
		Flux_Tensor[2][2] = Flux_Tensor[0][2]*val_velocity[1];
		Flux_Tensor[3][2] = Flux_Tensor[0][2]*val_velocity[2]+val_pressure; 
		Flux_Tensor[4][2] = Flux_Tensor[0][2]*val_enthalpy;

	}
	if(nDim == 2) {
		Flux_Tensor[0][0] = val_density*val_velocity[0]; 
		Flux_Tensor[1][0] = Flux_Tensor[0][0]*val_velocity[0]+val_pressure; 
		Flux_Tensor[2][0] = Flux_Tensor[0][0]*val_velocity[1];
		Flux_Tensor[3][0] = Flux_Tensor[0][0]*val_enthalpy;

		Flux_Tensor[0][1] = val_density*val_velocity[1];
		Flux_Tensor[1][1] = Flux_Tensor[0][1]*val_velocity[0];
		Flux_Tensor[2][1] = Flux_Tensor[0][1]*val_velocity[1]+val_pressure; 
		Flux_Tensor[3][1] = Flux_Tensor[0][1]*val_enthalpy;
	}
}

void CNumerics::GetInviscidProjFlux(double *val_density, double *val_velocity,
		double *val_pressure, double *val_enthalpy,
		double *val_normal, double *val_Proj_Flux) {
	if (nDim == 3) {
		double rhou = *val_density*val_velocity[0]; 
		double rhov = *val_density*val_velocity[1];
		double rhow = *val_density*val_velocity[2];

		val_Proj_Flux[0] = rhou*val_normal[0]; 
		val_Proj_Flux[1] = (rhou*val_velocity[0]+*val_pressure)*val_normal[0]; 
		val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
		val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
		val_Proj_Flux[4] = rhou**val_enthalpy*val_normal[0];

		val_Proj_Flux[0] += rhov*val_normal[1];
		val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
		val_Proj_Flux[2] += (rhov*val_velocity[1]+*val_pressure)*val_normal[1]; 
		val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
		val_Proj_Flux[4] += rhov**val_enthalpy*val_normal[1];

		val_Proj_Flux[0] += rhow*val_normal[2];
		val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
		val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
		val_Proj_Flux[3] += (rhow*val_velocity[2]+*val_pressure)*val_normal[2]; 
		val_Proj_Flux[4] += rhow**val_enthalpy*val_normal[2];
	}

	if(nDim == 2) {
		double rhou = *val_density*val_velocity[0]; 
		double rhov = *val_density*val_velocity[1];
		val_Proj_Flux[0] = rhou*val_normal[0] + rhov*val_normal[1]; 
		val_Proj_Flux[1] = (rhou*val_velocity[0]+*val_pressure)*val_normal[0] + rhov*val_velocity[0]*val_normal[1]; 
		val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0] + (rhov*val_velocity[1]+*val_pressure)*val_normal[1]; 
		val_Proj_Flux[3] = rhou**val_enthalpy*val_normal[0] + rhov**val_enthalpy*val_normal[1]; 
	}
}



void CNumerics::GetInviscidArtCompProjFlux(double *val_density, double *val_velocity,
		double *val_pressure, double *val_gravityforce, double *val_betainc2, double *val_normal,
		double *val_Proj_Flux) {
	if (nDim == 3) {
		double u = val_velocity[0]; 
		double v = val_velocity[1];
		double w = val_velocity[2];
		double rho = *val_density;
		double press = *val_pressure;
		double gravityforce = *val_gravityforce;
		double betainc2 = *val_betainc2;

		val_Proj_Flux[0] = betainc2*u*val_normal[0]; 
		val_Proj_Flux[1] = (u*u+(press/rho)+gravityforce)*val_normal[0]; 
		val_Proj_Flux[2] = u*v*val_normal[0];
		val_Proj_Flux[3] = u*w*val_normal[0];

		val_Proj_Flux[0] += betainc2*v*val_normal[1];
		val_Proj_Flux[1] += v*u*val_normal[1];
		val_Proj_Flux[2] += (v*v+(press/rho)+gravityforce)*val_normal[1]; 
		val_Proj_Flux[3] += v*w*val_normal[1];

		val_Proj_Flux[0] += betainc2*w*val_normal[2];
		val_Proj_Flux[1] += w*u*val_normal[2];
		val_Proj_Flux[2] += w*v*val_normal[2];
		val_Proj_Flux[3] += (w*w+(press/rho)+gravityforce)*val_normal[2]; 
	}

	if(nDim == 2) {
		double u = val_velocity[0]; 
		double v = val_velocity[1];
		double rho = *val_density;
		double gravityforce = *val_gravityforce;
		double press = *val_pressure;
		double betainc2 = *val_betainc2;

		val_Proj_Flux[0] = betainc2*u*val_normal[0] + betainc2*v*val_normal[1]; 
		val_Proj_Flux[1] = (u*u+(press/rho)+gravityforce)*val_normal[0] + v*u*val_normal[1]; 
		val_Proj_Flux[2] = u*v*val_normal[0] + (v*v+(press/rho)+gravityforce)*val_normal[1]; 

	}

}

void CNumerics::GetInviscidProjFlux(double *val_density, double **val_velocity,
		double *val_pressure, double *val_enthalpy,
		double *val_normal, double *val_Proj_Flux) {
	unsigned short iFluids, loc = 0;

	for ( iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = iFluids*(nDim+2);

		if(nDim == 3) {
			double rhou		= val_density[iFluids]*val_velocity[iFluids][0];
			double rhov		= val_density[iFluids]*val_velocity[iFluids][1];
			double rhow		= val_density[iFluids]*val_velocity[iFluids][2];
			double pressure = val_pressure[iFluids];
			double u		= val_velocity[iFluids][0];
			double v		= val_velocity[iFluids][1];
			double w		= val_velocity[iFluids][2];
			double enthalpy = val_enthalpy[iFluids];
			val_Proj_Flux[loc + 0] = rhou*val_normal[0] + rhov*val_normal[1] + rhow*val_normal[2];
			val_Proj_Flux[loc + 1] = (rhou*u + pressure)*val_normal[0] + rhov*u*val_normal[1]+ rhow*u*val_normal[1];
			val_Proj_Flux[loc + 2] = rhou*v*val_normal[0] + (rhov*v+pressure)*val_normal[1] + rhow*v*val_normal[0];
			val_Proj_Flux[loc + 3] = rhou*w*val_normal[0] + rhov*w*val_normal[1] + (rhow*w + pressure)*val_normal[2];
			val_Proj_Flux[loc + 4] = rhou*enthalpy*val_normal[0] + rhov*enthalpy*val_normal[1] + rhow*enthalpy*val_normal[2] ;
		}

		if(nDim == 2) {
			double rhou		= val_density[iFluids]*val_velocity[iFluids][0]; 
			double rhov		= val_density[iFluids]*val_velocity[iFluids][1];
			double pressure = val_pressure[iFluids];
			double u		= val_velocity[iFluids][0];
			double v		= val_velocity[iFluids][1];
			double enthalpy = val_enthalpy[iFluids];
			val_Proj_Flux[loc + 0] = rhou*val_normal[0] + rhov*val_normal[1]; 
			val_Proj_Flux[loc + 1] = (rhou*u + pressure)*val_normal[0] + rhov*u*val_normal[1]; 
			val_Proj_Flux[loc + 2] = rhou*v*val_normal[0] + (rhov*v+pressure)*val_normal[1]; 
			val_Proj_Flux[loc + 3] = rhou*enthalpy*val_normal[0] + rhov*enthalpy*val_normal[1]; 			
		}

	}
}

void CNumerics::GetInviscidProjFlux_(double *val_density, double **val_velocity,
		double *val_pressure, double *val_enthalpy,
		double *val_energy_vib, double *val_normal,
		double *val_Proj_Flux) {
	unsigned short iSpecies, loc = 0;

	
	
	for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	

		
		if (nDim == 3) {
			double rhou		= val_density[iSpecies]*val_velocity[iSpecies][0]; 
			double rhov		= val_density[iSpecies]*val_velocity[iSpecies][1];
			double rhow   = val_density[iSpecies]*val_velocity[iSpecies][2];
			double pressure = val_pressure[iSpecies];
			double u		= val_velocity[iSpecies][0];
			double v		= val_velocity[iSpecies][1];
			double w		= val_velocity[iSpecies][2];
			double enthalpy = val_enthalpy[iSpecies];
			double energy_vibu = val_energy_vib[iSpecies]*val_velocity[iSpecies][0];
			double energy_vibv = val_energy_vib[iSpecies]*val_velocity[iSpecies][1];
			double energy_vibw = val_energy_vib[iSpecies]*val_velocity[iSpecies][2];

			val_Proj_Flux[loc + 0] = rhou*val_normal[0]; 
			val_Proj_Flux[loc + 1] = (rhou*u + pressure)*val_normal[0];
			val_Proj_Flux[loc + 2] = rhou*v*val_normal[0];
			val_Proj_Flux[loc + 3] = rhou*w*val_normal[0];
			val_Proj_Flux[loc + 4] = rhou*enthalpy*val_normal[0];
			if (iSpecies < nDiatomics)
				val_Proj_Flux[loc + 5] = energy_vibu*val_normal[0];

			val_Proj_Flux[loc + 0] += rhov*val_normal[1]; 
			val_Proj_Flux[loc + 1] += rhov*u*val_normal[1];
			val_Proj_Flux[loc + 2] += (rhov*v + pressure)*val_normal[1];
			val_Proj_Flux[loc + 3] += rhov*w*val_normal[1];
			val_Proj_Flux[loc + 4] += rhov*enthalpy*val_normal[1];
			if (iSpecies < nDiatomics)
				val_Proj_Flux[loc + 5] += energy_vibv*val_normal[1];

			val_Proj_Flux[loc + 0] += rhow*val_normal[2]; 
			val_Proj_Flux[loc + 1] += rhow*u*val_normal[2];
			val_Proj_Flux[loc + 2] += rhow*v*val_normal[2];
			val_Proj_Flux[loc + 3] += (rhow*w + pressure)*val_normal[2];
			val_Proj_Flux[loc + 4] += rhow*enthalpy*val_normal[2];
			if (iSpecies < nDiatomics)
				val_Proj_Flux[loc + 5] += energy_vibw*val_normal[2];
		}
		
		if (nDim == 2) {
			
			
			double rhou		= val_density[iSpecies]*val_velocity[iSpecies][0]; 
			double rhov		= val_density[iSpecies]*val_velocity[iSpecies][1];
			double pressure = val_pressure[iSpecies];
			double u		= val_velocity[iSpecies][0];
			double v		= val_velocity[iSpecies][1];
			double enthalpy = val_enthalpy[iSpecies];
			double energy_vibu = val_energy_vib[iSpecies]*val_velocity[iSpecies][0];
			double energy_vibv = val_energy_vib[iSpecies]*val_velocity[iSpecies][1];
			
			val_Proj_Flux[loc + 0] = rhou*val_normal[0]; 
			val_Proj_Flux[loc + 1] = (rhou*u + pressure)*val_normal[0];
			val_Proj_Flux[loc + 2] = rhou*v*val_normal[0];
			val_Proj_Flux[loc + 3] = rhou*enthalpy*val_normal[0];
			if (iSpecies < nDiatomics)
				val_Proj_Flux[loc + 4] = energy_vibu*val_normal[0];
			
			val_Proj_Flux[loc + 0] += rhov*val_normal[1]; 
			val_Proj_Flux[loc + 1] += rhov*u*val_normal[1];
			val_Proj_Flux[loc + 2] += (rhov*v + pressure)*val_normal[1];
			val_Proj_Flux[loc + 3] += rhov*enthalpy*val_normal[1];
			if (iSpecies < nDiatomics)
				val_Proj_Flux[loc + 4] += energy_vibv*val_normal[1];						
			
		}
	}
}

void CNumerics::GetInviscidProjJac(double *val_velocity, double val_energy, double *val_normal, 
		double val_scale, double **val_Proj_Jac_Tensor) {
	unsigned short iDim, jDim;
	double sqvel = 0.0;
	double proj_vel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		sqvel    += val_velocity[iDim]*val_velocity[iDim];
		proj_vel += val_velocity[iDim]*val_normal[iDim];
	}
	double phi = 0.5*Gamma_Minus_One*sqvel;
	double a1 = Gamma*val_energy-phi;
	double a2 = Gamma-1.0;

	val_Proj_Jac_Tensor[0][0] = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		val_Proj_Jac_Tensor[0][iDim+1] = val_scale*val_normal[iDim];
	val_Proj_Jac_Tensor[0][nDim+1] = 0.0;

	for (iDim = 0; iDim < nDim; iDim++) {
		val_Proj_Jac_Tensor[iDim+1][0] = val_scale*(val_normal[iDim]*phi - val_velocity[iDim]*proj_vel);
		for (jDim = 0; jDim < nDim; jDim++)
			val_Proj_Jac_Tensor[iDim+1][jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iDim]-a2*val_normal[iDim]*val_velocity[jDim]);
		val_Proj_Jac_Tensor[iDim+1][iDim+1] += val_scale*proj_vel;
		val_Proj_Jac_Tensor[iDim+1][nDim+1] = val_scale*a2*val_normal[iDim];
	}

	val_Proj_Jac_Tensor[nDim+1][0] = val_scale*proj_vel*(phi-a1);
	for (iDim = 0; iDim < nDim; iDim++)
		val_Proj_Jac_Tensor[nDim+1][iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iDim]*proj_vel);
	val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*Gamma*proj_vel;
}

void CNumerics::GetInviscidArtCompProjJac(double val_density, double *val_velocity, double val_betainc2, double *val_normal, 
		double val_scale, double **val_Proj_Jac_Tensor) {
	unsigned short iDim;
	double proj_vel, sqvel;

	if (nDim == 2) {
		sqvel = 0.0; proj_vel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			proj_vel += val_velocity[iDim]*val_normal[iDim];

		val_Proj_Jac_Tensor[0][0] = 0.0;
		val_Proj_Jac_Tensor[0][1] = val_scale*val_betainc2*val_normal[0];
		val_Proj_Jac_Tensor[0][2] = val_scale*val_betainc2*val_normal[1];

		val_Proj_Jac_Tensor[1][0] = val_scale*val_normal[0]/val_density;
		val_Proj_Jac_Tensor[1][1] = val_scale*(val_velocity[0]*val_normal[0] + proj_vel);
		val_Proj_Jac_Tensor[1][2] = val_scale*val_velocity[0]*val_normal[1];

		val_Proj_Jac_Tensor[2][0] = val_scale*val_normal[1]/val_density;
		val_Proj_Jac_Tensor[2][1] = val_scale*val_velocity[1]*val_normal[0];
		val_Proj_Jac_Tensor[2][2] = val_scale*(val_velocity[1]*val_normal[1] + proj_vel);
	}

}

void CNumerics::GetInviscidProjJac(double **val_velocity, double *val_energy, double *val_normal, 
		double val_scale, double **val_Proj_Jac_Tensor) {
	unsigned short iDim, jDim;
	unsigned short iVar, jVar, iFluids, loc = 0;

	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++) {
			val_Proj_Jac_Tensor[iVar][jVar] = 0.0;
		}
	}
	double sqvel,proj_vel;
	double phi, a1, a2;

	for(iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = iFluids * (nDim+2);
		sqvel = 0.0;
		proj_vel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			sqvel    += val_velocity[iFluids][iDim]*val_velocity[iFluids][iDim];
			proj_vel += val_velocity[iFluids][iDim]*val_normal[iDim];
		}
		phi = 0.5*Gamma_Minus_One*sqvel;
		a1  = Gamma*val_energy[iFluids]-phi;
		a2  = Gamma-1.0;

		val_Proj_Jac_Tensor[loc+0][loc+0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			val_Proj_Jac_Tensor[loc+0][loc+iDim+1] = val_scale*val_normal[iDim];
		val_Proj_Jac_Tensor[loc+0][loc+nDim+1] = 0.0;

		for (iDim = 0; iDim < nDim; iDim++) {
			val_Proj_Jac_Tensor[loc+iDim+1][loc+0] = val_scale*(val_normal[iDim]*phi - val_velocity[iFluids][iDim]*proj_vel);
			for (jDim = 0; jDim < nDim; jDim++)
				val_Proj_Jac_Tensor[loc+iDim+1][loc+jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iFluids][iDim]-a2*val_normal[iDim]*val_velocity[iFluids][jDim]);
			val_Proj_Jac_Tensor[loc+iDim+1][loc+iDim+1] += val_scale*proj_vel;
			val_Proj_Jac_Tensor[loc+iDim+1][loc+nDim+1] = val_scale*a2*val_normal[iDim];
		}

		val_Proj_Jac_Tensor[loc+nDim+1][loc+0] = val_scale*proj_vel*(phi-a1);
		for (iDim = 0; iDim < nDim; iDim++)
			val_Proj_Jac_Tensor[loc+nDim+1][loc+iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iFluids][iDim]*proj_vel);
		val_Proj_Jac_Tensor[loc+nDim+1][loc+nDim+1] = val_scale*Gamma*proj_vel;
	}
}

void CNumerics::GetInviscidProjJac_(double **val_velocity, double *val_energy, double *val_normal,
		double val_scale, double **val_Proj_Jac_Tensor) {
	unsigned short iDim, jDim;
	unsigned short iVar, jVar, iSpecies, loc = 0;

	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++) {
			val_Proj_Jac_Tensor[iVar][jVar] = 0.0;
		}
	}
	double sqvel,proj_vel;
	double phi, a1, a2;

	for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		sqvel = 0.0;
		proj_vel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			sqvel    += val_velocity[iSpecies][iDim]*val_velocity[iSpecies][iDim];
			proj_vel += val_velocity[iSpecies][iDim]*val_normal[iDim];
		}
		if (iSpecies < nDiatomics){
			phi = 0.5*(GammaDiatomic-1.0)*sqvel;
			a1  = GammaDiatomic*val_energy[iSpecies]-phi;
			a2  = GammaDiatomic-1.0;
		}
		else {
			phi = 0.5*(GammaMonatomic-1.0)*sqvel;
			a1  = GammaMonatomic*val_energy[iSpecies]-phi;
			a2  = GammaMonatomic-1.0;
		}

		val_Proj_Jac_Tensor[loc+0][loc+0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			val_Proj_Jac_Tensor[loc+0][loc+iDim+1] = val_scale*val_normal[iDim];
		val_Proj_Jac_Tensor[loc+0][loc+nDim+1] = 0.0;

		for (iDim = 0; iDim < nDim; iDim++) {
			val_Proj_Jac_Tensor[loc+iDim+1][loc+0] = val_scale*(val_normal[iDim]*phi - val_velocity[iSpecies][iDim]*proj_vel);
			for (jDim = 0; jDim < nDim; jDim++)
				val_Proj_Jac_Tensor[loc+iDim+1][loc+jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iSpecies][iDim]-a2*val_normal[iDim]*val_velocity[iSpecies][jDim]);
			val_Proj_Jac_Tensor[loc+iDim+1][loc+iDim+1] += val_scale*proj_vel;
			val_Proj_Jac_Tensor[loc+iDim+1][loc+nDim+1] = val_scale*a2*val_normal[iDim];
		}

		val_Proj_Jac_Tensor[loc+nDim+1][loc+0] = val_scale*proj_vel*(phi-a1);
		for (iDim = 0; iDim < nDim; iDim++)
			val_Proj_Jac_Tensor[loc+nDim+1][loc+iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iSpecies][iDim]*proj_vel);
		val_Proj_Jac_Tensor[loc+nDim+1][loc+nDim+1] = val_scale*Gamma*proj_vel;
	}
}


void CNumerics::SetPastSol (double *val_u_nM1, double *val_u_n, double *val_u_nP1) {
	unsigned short iVar;

	for(iVar = 0; iVar < nVar; iVar++) {
		U_nM1[iVar] = val_u_nM1[iVar];
		U_n[iVar] = val_u_n[iVar];
		U_nP1[iVar] = val_u_nP1[iVar];
	}

}
void CNumerics::SetPastVolume (double val_volume_nM1, double val_volume_n, double val_volume_nP1) {
	Volume_nM1 = val_volume_nM1;
	Volume_n = val_volume_n;
	Volume_nP1 = val_volume_nP1;
}
void CNumerics::SetTimeStep (double val_timestep) {
	TimeStep = val_timestep;
}

void CNumerics::GetPMatrix(double *val_density, double *val_velocity,
		double *val_soundspeed, double *val_normal, double **val_p_tensor) {
	double sqvel, rhooc, rhoxc, c2;

	rhooc = *val_density / *val_soundspeed, 
			rhoxc = *val_density * *val_soundspeed,
			c2 = *val_soundspeed * *val_soundspeed;

	if (nDim == 3) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

		val_p_tensor[0][0]=val_normal[0];
		val_p_tensor[0][1]=val_normal[1];
		val_p_tensor[0][2]=val_normal[2];
		val_p_tensor[0][3]=0.5*rhooc;
		val_p_tensor[0][4]=0.5*rhooc;	

		val_p_tensor[1][0]=val_velocity[0]*val_normal[0];
		val_p_tensor[1][1]=val_velocity[0]*val_normal[1]-*val_density*val_normal[2];
		val_p_tensor[1][2]=val_velocity[0]*val_normal[2]+*val_density*val_normal[1];
		val_p_tensor[1][3]=0.5*(val_velocity[0]*rhooc+*val_density*val_normal[0]);
		val_p_tensor[1][4]=0.5*(val_velocity[0]*rhooc-*val_density*val_normal[0]);	

		val_p_tensor[2][0]=val_velocity[1]*val_normal[0]+*val_density*val_normal[2];
		val_p_tensor[2][1]=val_velocity[1]*val_normal[1];
		val_p_tensor[2][2]=val_velocity[1]*val_normal[2]-*val_density*val_normal[0];
		val_p_tensor[2][3]=0.5*(val_velocity[1]*rhooc+*val_density*val_normal[1]);
		val_p_tensor[2][4]=0.5*(val_velocity[1]*rhooc-*val_density*val_normal[1]);

		val_p_tensor[3][0]=val_velocity[2]*val_normal[0]-*val_density*val_normal[1];
		val_p_tensor[3][1]=val_velocity[2]*val_normal[1]+*val_density*val_normal[0];
		val_p_tensor[3][2]=val_velocity[2]*val_normal[2];
		val_p_tensor[3][3]=0.5*(val_velocity[2]*rhooc+*val_density*val_normal[2]);
		val_p_tensor[3][4]=0.5*(val_velocity[2]*rhooc-*val_density*val_normal[2]);	

		val_p_tensor[4][0]=0.5*sqvel*val_normal[0]+*val_density*val_velocity[1]*val_normal[2]-*val_density*val_velocity[2]*val_normal[1];
		val_p_tensor[4][1]=0.5*sqvel*val_normal[1]-*val_density*val_velocity[0]*val_normal[2]+*val_density*val_velocity[2]*val_normal[0];
		val_p_tensor[4][2]=0.5*sqvel*val_normal[2]+*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
		val_p_tensor[4][3]=0.5*(0.5*sqvel*rhooc+*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2])+rhoxc/Gamma_Minus_One);
		val_p_tensor[4][4]=0.5*(0.5*sqvel*rhooc-*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2])+rhoxc/Gamma_Minus_One);
	}

	if(nDim == 2) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

		val_p_tensor[0][0]=1.0;
		val_p_tensor[0][1]=0.0;
		val_p_tensor[0][2]=0.5*rhooc;
		val_p_tensor[0][3]=0.5*rhooc;	

		val_p_tensor[1][0]=val_velocity[0];
		val_p_tensor[1][1]=*val_density*val_normal[1];
		val_p_tensor[1][2]=0.5*(val_velocity[0]*rhooc+val_normal[0]**val_density);
		val_p_tensor[1][3]=0.5*(val_velocity[0]*rhooc-val_normal[0]**val_density);	

		val_p_tensor[2][0]=val_velocity[1];
		val_p_tensor[2][1]=-*val_density*val_normal[0];
		val_p_tensor[2][2]=0.5*(val_velocity[1]*rhooc+val_normal[1]**val_density);
		val_p_tensor[2][3]=0.5*(val_velocity[1]*rhooc-val_normal[1]**val_density);

		val_p_tensor[3][0]=0.5*sqvel;
		val_p_tensor[3][1]=*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
		val_p_tensor[3][2]=0.5*(0.5*sqvel*rhooc+*val_density*val_velocity[0]*val_normal[0]+*val_density*val_velocity[1]*val_normal[1]+rhoxc/Gamma_Minus_One);
		val_p_tensor[3][3]=0.5*(0.5*sqvel*rhooc-*val_density*val_velocity[0]*val_normal[0]-*val_density*val_velocity[1]*val_normal[1]+rhoxc/Gamma_Minus_One);
	}
}

void CNumerics::GetPMatrix(double *val_density, double **val_velocity,
		double *val_soundspeed, double *val_normal, double **val_p_tensor) {

	double sqvel, rhooc, rhoxc, c2, u ,v, w,rho ;
	unsigned short loc, iVar, jVar, iFluids;

	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++) {
			val_p_tensor[iVar][jVar]= 0.0;
		}
	}

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		rhooc = val_density[iFluids]    / val_soundspeed[iFluids];
		rhoxc = val_density[iFluids]	* val_soundspeed[iFluids];
		c2    = val_soundspeed[iFluids] * val_soundspeed[iFluids];
		rho   = val_density[iFluids];
		loc   = iFluids * (nDim+2);

		if (nDim == 3) {
			u	  = val_velocity[iFluids][0];
			v     = val_velocity[iFluids][1];
			w     = val_velocity[iFluids][2];
			sqvel = u*u + v*v + w*w;

			val_p_tensor[loc +0][loc +0]=val_normal[0];
			val_p_tensor[loc +0][loc +1]=val_normal[1];
			val_p_tensor[loc +0][loc +2]=val_normal[2];
			val_p_tensor[loc +0][loc +3]=0.5*rhooc;
			val_p_tensor[loc +0][loc +4]=0.5*rhooc;

			val_p_tensor[loc +1][loc +0]=u*val_normal[0];
			val_p_tensor[loc +1][loc +1]=u*val_normal[1]-rho*val_normal[2];
			val_p_tensor[loc +1][loc +2]=u*val_normal[2]+rho*val_normal[1];
			val_p_tensor[loc +1][loc +3]=0.5*(u*rhooc+rho*val_normal[0]);
			val_p_tensor[loc +1][loc +4]=0.5*(u*rhooc-rho*val_normal[0]);

			val_p_tensor[loc +2][loc +0]=v*val_normal[0]+rho*val_normal[2];
			val_p_tensor[loc +2][loc +1]=v*val_normal[1];
			val_p_tensor[loc +2][loc +2]=v*val_normal[2]-rho*val_normal[0];
			val_p_tensor[loc +2][loc +3]=0.5*(v*rhooc+rho*val_normal[1]);
			val_p_tensor[loc +2][loc +4]=0.5*(v*rhooc-rho*val_normal[1]);

			val_p_tensor[loc +3][loc +0]=w*val_normal[0]-rho*val_normal[1];
			val_p_tensor[loc +3][loc +1]=w*val_normal[1]+rho*val_normal[0];
			val_p_tensor[loc +3][loc +2]=w*val_normal[2];
			val_p_tensor[loc +3][loc +3]=0.5*(w*rhooc+rho*val_normal[2]);
			val_p_tensor[loc +3][loc +4]=0.5*(w*rhooc-rho*val_normal[2]);

			val_p_tensor[loc +4][loc +0]=0.5*sqvel*val_normal[0]+rho*v*val_normal[2]-rho*w*val_normal[1];
			val_p_tensor[loc +4][loc +1]=0.5*sqvel*val_normal[1]-rho*u*val_normal[2]+rho*w*val_normal[0];
			val_p_tensor[loc +4][loc +2]=0.5*sqvel*val_normal[2]+rho*u*val_normal[1]-rho*v*val_normal[0];
			val_p_tensor[loc +4][loc +3]=0.5*(0.5*sqvel*rhooc+rho*(u*val_normal[0]+v*val_normal[1]+w*val_normal[2])+rhoxc/Gamma_Minus_One);
			val_p_tensor[loc +4][loc +4]=0.5*(0.5*sqvel*rhooc-rho*(u*val_normal[0]+v*val_normal[1]+w*val_normal[2])+rhoxc/Gamma_Minus_One);
		}
		if(nDim == 2) {
			u	  = val_velocity[iFluids][0];
			v     = val_velocity[iFluids][1];
			sqvel = u*u + v*v;

			val_p_tensor[loc+0][loc+0]=1.0;
			val_p_tensor[loc+0][loc+1]=0.0;
			val_p_tensor[loc+0][loc+2]=0.5*rhooc;
			val_p_tensor[loc+0][loc+3]=0.5*rhooc;	

			val_p_tensor[loc+1][loc+0]=u;
			val_p_tensor[loc+1][loc+1]=rho *val_normal[1];
			val_p_tensor[loc+1][loc+2]=0.5*(u*rhooc+val_normal[0]*rho);
			val_p_tensor[loc+1][loc+3]=0.5*(u*rhooc-val_normal[0]*rho);	

			val_p_tensor[loc+2][loc+0]=v;
			val_p_tensor[loc+2][loc+1]=-rho*val_normal[0];
			val_p_tensor[loc+2][loc+2]=0.5*(v*rhooc+val_normal[1]*rho);
			val_p_tensor[loc+2][loc+3]=0.5*(v*rhooc-val_normal[1]*rho);

			val_p_tensor[loc+3][loc+0]=0.5*sqvel;
			val_p_tensor[loc+3][loc+1]=rho*u*val_normal[1]-rho*v*val_normal[0];
			val_p_tensor[loc+3][loc+2]=0.5*(0.5*sqvel*rhooc+rho*u*val_normal[0]+rho*v*val_normal[1]+rhoxc/Gamma_Minus_One);
			val_p_tensor[loc+3][loc+3]=0.5*(0.5*sqvel*rhooc-rho*u*val_normal[0]-rho*v*val_normal[1]+rhoxc/Gamma_Minus_One);

		}
	}
}

void CNumerics::GetPMatrix_(double *val_density, double **val_velocity,
														double *val_soundspeed, double *val_normal, double **val_p_tensor) {
	
	double sqvel, rhooc, rhoxc, c2;
	unsigned short loc, iVar, jVar, iSpecies;
	
	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++) {
			val_p_tensor[iVar][jVar]= 0.0;
		}
	}
	
	if (nDim == 3) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {			
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	
			
			rhooc = val_density[iSpecies]	 / val_soundspeed[iSpecies], 
			rhoxc = val_density[iSpecies]	 * val_soundspeed[iSpecies], 
			c2 =    val_soundspeed[iSpecies] * val_soundspeed[iSpecies];				
			sqvel = val_velocity[iSpecies][0]*val_velocity[iSpecies][0]+val_velocity[iSpecies][1]*val_velocity[iSpecies][1]+val_velocity[iSpecies][2]*val_velocity[iSpecies][2];
			
			val_p_tensor[loc+0][loc+0]=val_normal[0];
			val_p_tensor[loc+0][loc+1]=val_normal[1];
			val_p_tensor[loc+0][loc+2]=val_normal[2];
			val_p_tensor[loc+0][loc+3]=0.5*rhooc;
			val_p_tensor[loc+0][loc+4]=0.5*rhooc;	
			
			val_p_tensor[loc+1][loc+0]=val_velocity[iSpecies][0]*val_normal[0];
			val_p_tensor[loc+1][loc+1]=val_velocity[iSpecies][0]*val_normal[1]-val_density[iSpecies]*val_normal[2];
			val_p_tensor[loc+1][loc+2]=val_velocity[iSpecies][0]*val_normal[2]+val_density[iSpecies]*val_normal[1];
			val_p_tensor[loc+1][loc+3]=0.5*(val_velocity[iSpecies][0]*rhooc+val_density[iSpecies]*val_normal[0]);
			val_p_tensor[loc+1][loc+4]=0.5*(val_velocity[iSpecies][0]*rhooc-val_density[iSpecies]*val_normal[0]);	
			
			val_p_tensor[loc+2][loc+0]=val_velocity[iSpecies][1]*val_normal[0]+val_density[iSpecies]*val_normal[2];
			val_p_tensor[loc+2][loc+1]=val_velocity[iSpecies][1]*val_normal[1];
			val_p_tensor[loc+2][loc+2]=val_velocity[iSpecies][1]*val_normal[2]-val_density[iSpecies]*val_normal[0];
			val_p_tensor[loc+2][loc+3]=0.5*(val_velocity[iSpecies][1]*rhooc+val_density[iSpecies]*val_normal[1]);
			val_p_tensor[loc+2][loc+4]=0.5*(val_velocity[iSpecies][1]*rhooc-val_density[iSpecies]*val_normal[1]);
			
			val_p_tensor[loc+3][loc+0]=val_velocity[iSpecies][2]*val_normal[0]-val_density[iSpecies]*val_normal[1];
			val_p_tensor[loc+3][loc+1]=val_velocity[iSpecies][2]*val_normal[1]+val_density[iSpecies]*val_normal[0];
			val_p_tensor[loc+3][loc+2]=val_velocity[iSpecies][2]*val_normal[2];
			val_p_tensor[loc+3][loc+3]=0.5*(val_velocity[iSpecies][2]*rhooc+val_density[iSpecies]*val_normal[2]);
			val_p_tensor[loc+3][loc+4]=0.5*(val_velocity[iSpecies][2]*rhooc-val_density[iSpecies]*val_normal[2]);	
			
			val_p_tensor[loc+4][loc+0]=0.5*sqvel*val_normal[0]+val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[2]-val_density[iSpecies]*val_velocity[iSpecies][2]*val_normal[1];
			val_p_tensor[loc+4][loc+1]=0.5*sqvel*val_normal[1]-val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[2]+val_density[iSpecies]*val_velocity[iSpecies][2]*val_normal[0];
			val_p_tensor[loc+4][loc+2]=0.5*sqvel*val_normal[2]+val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[1]-val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[0];
			if (iSpecies < nDiatomics) {
				val_p_tensor[loc+4][loc+3]=0.5*(0.5*sqvel*rhooc+val_density[iSpecies]*(val_velocity[iSpecies][0]*val_normal[0]+val_velocity[iSpecies][1]*val_normal[1]+val_velocity[iSpecies][2]*val_normal[2])+rhoxc/(GammaDiatomic-1.0));
				val_p_tensor[loc+4][loc+4]=0.5*(0.5*sqvel*rhooc-val_density[iSpecies]*(val_velocity[iSpecies][0]*val_normal[0]+val_velocity[iSpecies][1]*val_normal[1]+val_velocity[iSpecies][2]*val_normal[2])+rhoxc/(GammaDiatomic-1.0));
			}
			else {
				val_p_tensor[loc+4][loc+3]=0.5*(0.5*sqvel*rhooc+val_density[iSpecies]*(val_velocity[iSpecies][0]*val_normal[0]+val_velocity[iSpecies][1]*val_normal[1]+val_velocity[iSpecies][2]*val_normal[2])+rhoxc/(GammaMonatomic-1.0));
				val_p_tensor[loc+4][loc+4]=0.5*(0.5*sqvel*rhooc-val_density[iSpecies]*(val_velocity[iSpecies][0]*val_normal[0]+val_velocity[iSpecies][1]*val_normal[1]+val_velocity[iSpecies][2]*val_normal[2])+rhoxc/(GammaMonatomic-1.0));
			}
			
		}
	}
	
	if(nDim == 2) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	
			
			rhooc = val_density[iSpecies]	 / val_soundspeed[iSpecies], 
			rhoxc = val_density[iSpecies]	 * val_soundspeed[iSpecies], 
			c2 =    val_soundspeed[iSpecies] * val_soundspeed[iSpecies];
			sqvel = val_velocity[iSpecies][0]*val_velocity[iSpecies][0]+val_velocity[iSpecies][1]*val_velocity[iSpecies][1];
			
			val_p_tensor[loc+0][loc+0]=1.0;
			val_p_tensor[loc+0][loc+1]=0.0;
			val_p_tensor[loc+0][loc+2]=0.5*rhooc;
			val_p_tensor[loc+0][loc+3]=0.5*rhooc;	
			
			val_p_tensor[loc+1][loc+0]=val_velocity[iSpecies][0];
			val_p_tensor[loc+1][loc+1]=val_density[iSpecies]*val_normal[1];
			val_p_tensor[loc+1][loc+2]=0.5*(val_velocity[iSpecies][0]*rhooc+val_normal[0]*val_density[iSpecies]);
			val_p_tensor[loc+1][loc+3]=0.5*(val_velocity[iSpecies][0]*rhooc-val_normal[0]*val_density[iSpecies]);	
			
			val_p_tensor[loc+2][loc+0]=val_velocity[iSpecies][1];
			val_p_tensor[loc+2][loc+1]=-val_density[iSpecies]*val_normal[0];
			val_p_tensor[loc+2][loc+2]=0.5*(val_velocity[iSpecies][1]*rhooc+val_normal[1]*val_density[iSpecies]);
			val_p_tensor[loc+2][loc+3]=0.5*(val_velocity[iSpecies][1]*rhooc-val_normal[1]*val_density[iSpecies]);
			
			val_p_tensor[loc+3][loc+0]=0.5*sqvel;
			val_p_tensor[loc+3][loc+1]=val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[1]-val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[0];
			if (iSpecies < nDiatomics) {
				val_p_tensor[loc+3][loc+2]=0.5*(0.5*sqvel*rhooc+val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[0]+val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[1]+rhoxc/(GammaDiatomic-1.0));
				val_p_tensor[loc+3][loc+3]=0.5*(0.5*sqvel*rhooc-val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[0]-val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[1]+rhoxc/(GammaDiatomic-1.0));
			}
			else {
				val_p_tensor[loc+3][loc+2]=0.5*(0.5*sqvel*rhooc+val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[0]+val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[1]+rhoxc/(GammaMonatomic-1.0));
				val_p_tensor[loc+3][loc+3]=0.5*(0.5*sqvel*rhooc-val_density[iSpecies]*val_velocity[iSpecies][0]*val_normal[0]-val_density[iSpecies]*val_velocity[iSpecies][1]*val_normal[1]+rhoxc/(GammaMonatomic-1.0));
			}
		}
	}
}

void CNumerics::GetPMatrix_inv(double *val_density, double *val_velocity,
		double *val_soundspeed, double *val_normal, double **val_invp_tensor) {
	double rhoxc, c2, gm1, k0orho, k1orho, gm1_o_c2, gm1_o_rhoxc, sqvel;

	rhoxc = *val_density * *val_soundspeed;
	c2 = *val_soundspeed * *val_soundspeed;
	gm1 = Gamma_Minus_One;
	k0orho = val_normal[0] / *val_density;
	k1orho = val_normal[1] / *val_density;
	gm1_o_c2 = gm1/c2;
	gm1_o_rhoxc = gm1/rhoxc;

	if (nDim == 3) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

		val_invp_tensor[0][0]=val_normal[0]-val_normal[2]*val_velocity[1] / *val_density+val_normal[1]*val_velocity[2] / *val_density-val_normal[0]*0.5*gm1*sqvel/c2;
		val_invp_tensor[0][1]=val_normal[0]*gm1*val_velocity[0]/c2;
		val_invp_tensor[0][2]=val_normal[2] / *val_density+val_normal[0]*gm1*val_velocity[1]/c2;
		val_invp_tensor[0][3]=-val_normal[1] / *val_density+val_normal[0]*gm1*val_velocity[2]/c2;
		val_invp_tensor[0][4]=-val_normal[0]*gm1/c2;

		val_invp_tensor[1][0]=val_normal[1]+val_normal[2]*val_velocity[0] / *val_density-val_normal[0]*val_velocity[2] / *val_density-val_normal[1]*0.5*gm1*sqvel/c2;
		val_invp_tensor[1][1]=-val_normal[2] / *val_density+val_normal[1]*gm1*val_velocity[0]/c2;
		val_invp_tensor[1][2]=val_normal[1]*gm1*val_velocity[1]/c2;
		val_invp_tensor[1][3]=val_normal[0] / *val_density+val_normal[1]*gm1*val_velocity[2]/c2;
		val_invp_tensor[1][4]=-val_normal[1]*gm1/c2;

		val_invp_tensor[2][0]=val_normal[2]-val_normal[1]*val_velocity[0] / *val_density+val_normal[0]*val_velocity[1] / *val_density-val_normal[2]*0.5*gm1*sqvel/c2;
		val_invp_tensor[2][1]=val_normal[1] / *val_density+val_normal[2]*gm1*val_velocity[0]/c2;
		val_invp_tensor[2][2]=-val_normal[0] / *val_density+val_normal[2]*gm1*val_velocity[1]/c2;
		val_invp_tensor[2][3]=val_normal[2]*gm1*val_velocity[2]/c2;
		val_invp_tensor[2][4]=-val_normal[2]*gm1/c2;

		val_invp_tensor[3][0]=-(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+0.5*gm1*sqvel/rhoxc;
		val_invp_tensor[3][1]=val_normal[0] / *val_density-gm1*val_velocity[0]/rhoxc;
		val_invp_tensor[3][2]=val_normal[1] / *val_density-gm1*val_velocity[1]/rhoxc;
		val_invp_tensor[3][3]=val_normal[2] / *val_density-gm1*val_velocity[2]/rhoxc;
		val_invp_tensor[3][4]=Gamma_Minus_One/rhoxc;

		val_invp_tensor[4][0]=(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+0.5*gm1*sqvel/rhoxc;
		val_invp_tensor[4][1]=-val_normal[0] / *val_density-gm1*val_velocity[0]/rhoxc;
		val_invp_tensor[4][2]=-val_normal[1] / *val_density-gm1*val_velocity[1]/rhoxc;
		val_invp_tensor[4][3]=-val_normal[2] / *val_density-gm1*val_velocity[2]/rhoxc;
		val_invp_tensor[4][4]=Gamma_Minus_One/rhoxc;
	}
	if(nDim == 2) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

		val_invp_tensor[0][0]=1.0-0.5*gm1_o_c2*sqvel;
		val_invp_tensor[0][1]=gm1_o_c2*val_velocity[0];
		val_invp_tensor[0][2]=gm1_o_c2*val_velocity[1];
		val_invp_tensor[0][3]=-gm1_o_c2;

		val_invp_tensor[1][0]=-k1orho*val_velocity[0]+k0orho*val_velocity[1];
		val_invp_tensor[1][1]=k1orho;
		val_invp_tensor[1][2]=-k0orho;
		val_invp_tensor[1][3]=0.0;

		val_invp_tensor[2][0]=-k0orho*val_velocity[0]-k1orho*val_velocity[1]+0.5*gm1_o_rhoxc*sqvel;
		val_invp_tensor[2][1]=k0orho-gm1_o_rhoxc*val_velocity[0];
		val_invp_tensor[2][2]=k1orho-gm1_o_rhoxc*val_velocity[1];
		val_invp_tensor[2][3]=gm1_o_rhoxc;

		val_invp_tensor[3][0]=k0orho*val_velocity[0]+k1orho*val_velocity[1]+0.5*gm1_o_rhoxc*sqvel;
		val_invp_tensor[3][1]=-k0orho-gm1_o_rhoxc*val_velocity[0];
		val_invp_tensor[3][2]=-k1orho-gm1_o_rhoxc*val_velocity[1];
		val_invp_tensor[3][3]=gm1_o_rhoxc;
	}
}

void CNumerics::GetPMatrix_inv(double *val_density, double **val_velocity,
		double *val_soundspeed, double *val_normal, double **val_invp_tensor) {

	double rhoxc, c2, gm1, k0orho, k1orho, gm1_o_c2, gm1_o_rhoxc, sqvel,u ,v,w, rho;
	double k0, k1 ,k2;
	unsigned short loc, iFluids; // location along the matrix
	unsigned short iVar, jVar;
	gm1 = Gamma_Minus_One;

	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++) {
			val_invp_tensor[iVar][jVar] = 0.0;
		}
	}

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc			= (nDim+2)*iFluids;
		rhoxc		= val_density[iFluids] * val_soundspeed[iFluids];
		c2			= val_soundspeed[iFluids] * val_soundspeed[iFluids];
		k0orho		= val_normal[0] / val_density[iFluids];
		k1orho		= val_normal[1] / val_density[iFluids];
		gm1_o_c2	= gm1/c2;
		gm1_o_rhoxc = gm1/rhoxc;
		rho 		= val_density[iFluids];

		if (nDim == 3) {
			u = val_velocity[iFluids][0];
			v = val_velocity[iFluids][1];
			w = val_velocity[iFluids][2];
			k0 = val_normal[0];
			k1 = val_normal[1];
			k2 = val_normal[2];
			sqvel 	= u*u + v*v + w*w;

			val_invp_tensor[loc+ 0][loc + 0] =  k0 - k2*v/rho + k1*w/rho - k0*0.5*gm1*sqvel/c2;
			val_invp_tensor[loc + 0][loc + 1] =  k0*gm1*u/c2;
			val_invp_tensor[loc + 0][loc + 2] =  k2/rho + k0*gm1*v/c2;
			val_invp_tensor[loc + 0][loc + 3] = -k1/rho + k0*gm1*w/c2;
			val_invp_tensor[loc + 0][loc + 4] = -k0*gm1/c2;

			val_invp_tensor[loc +1][loc + 0] =  k1+k2*u/rho - k0*w/rho - k1*0.5*gm1*sqvel/c2;
			val_invp_tensor[loc +1][loc + 1] = -k2/rho + k1*gm1*u/c2;
			val_invp_tensor[loc +1][loc + 2] =  k1*gm1*v/c2;
			val_invp_tensor[loc +1][loc + 3] =  k0/rho + k1*gm1*w/c2;
			val_invp_tensor[loc +1][loc + 4] = -k1*gm1/c2;

			val_invp_tensor[loc +2][loc + 0] =  k2 - k1*u/rho + k0*v/rho - k2*0.5*gm1*sqvel/c2;
			val_invp_tensor[loc +2][loc + 1] =  k1/rho + k2*gm1*u/c2;
			val_invp_tensor[loc +2][loc + 2] = -k0/rho + k2*gm1*v/c2;
			val_invp_tensor[loc +2][loc + 3] =  k2*gm1*w/c2;
			val_invp_tensor[loc +2][loc + 4] = -k2*gm1/c2;

			val_invp_tensor[loc + 3][loc + 0] = -(k0*u+k1*v+k2*w)/rho + 0.5*gm1*sqvel/rhoxc;
			val_invp_tensor[loc + 3][loc + 1] = k0/rho - gm1*u/rhoxc;
			val_invp_tensor[loc + 3][loc + 2] = k1/rho - gm1*v/rhoxc;
			val_invp_tensor[loc + 3][loc + 3] = k2/rho - gm1*w/rhoxc;
			val_invp_tensor[loc + 3][loc + 4] = gm1/rhoxc;

			val_invp_tensor[loc + 4][loc + 0] = (k0*u+k1*v+k2*w) /rho + 0.5*gm1*sqvel/rhoxc;
			val_invp_tensor[loc + 4][loc + 1] = -k0/rho - gm1*u/rhoxc;
			val_invp_tensor[loc + 4][loc + 2] = -k1/rho - gm1*v/rhoxc;
			val_invp_tensor[loc + 4][loc + 3] = -k2/rho - gm1*w/rhoxc;
			val_invp_tensor[loc + 4][loc + 4] =  gm1/rhoxc;
		}

		if(nDim == 2) {

			sqvel = val_velocity[iFluids][0]*val_velocity[iFluids][0]+val_velocity[iFluids][1]*val_velocity[iFluids][1];

			val_invp_tensor[loc+0][loc+0]=1.0-0.5*gm1_o_c2*sqvel;
			val_invp_tensor[loc+0][loc+1]=gm1_o_c2*val_velocity[iFluids][0];
			val_invp_tensor[loc+0][loc+2]=gm1_o_c2*val_velocity[iFluids][1];
			val_invp_tensor[loc+0][loc+3]=-gm1_o_c2;

			val_invp_tensor[loc+1][loc+0]=-k1orho*val_velocity[iFluids][0]+k0orho*val_velocity[iFluids][1];
			val_invp_tensor[loc+1][loc+1]=k1orho;
			val_invp_tensor[loc+1][loc+2]=-k0orho;
			val_invp_tensor[loc+1][loc+3]=0.0;

			val_invp_tensor[loc+2][loc+0]=-k0orho*val_velocity[iFluids][0]-k1orho*val_velocity[iFluids][1]+0.5*gm1_o_rhoxc*sqvel;
			val_invp_tensor[loc+2][loc+1]=k0orho-gm1_o_rhoxc*val_velocity[iFluids][0];
			val_invp_tensor[loc+2][loc+2]=k1orho-gm1_o_rhoxc*val_velocity[iFluids][1];
			val_invp_tensor[loc+2][loc+3]=gm1_o_rhoxc;

			val_invp_tensor[loc+3][loc+0]=k0orho*val_velocity[iFluids][0]+k1orho*val_velocity[iFluids][1]+0.5*gm1_o_rhoxc*sqvel;
			val_invp_tensor[loc+3][loc+1]=-k0orho-gm1_o_rhoxc*val_velocity[iFluids][0];
			val_invp_tensor[loc+3][loc+2]=-k1orho-gm1_o_rhoxc*val_velocity[iFluids][1];
			val_invp_tensor[loc+3][loc+3]=gm1_o_rhoxc;
		}
	}
}

void CNumerics::GetPMatrix_inv_(double *val_density, double **val_velocity,
		double *val_soundspeed, double *val_normal, double **val_invp_tensor) {

	double rhoxc, c2, gm1, k0orho, k1orho, gm1_o_c2, gm1_o_rhoxc, sqvel;
	unsigned short loc, iSpecies; // location along the matrix
	unsigned short iVar, jVar;

	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++) {
			val_invp_tensor[iVar][jVar] = 0.0;
		}
	}

	if (nDim == 3) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

			rhoxc		= val_density[iSpecies] * val_soundspeed[iSpecies];
			c2			= val_soundspeed[iSpecies] * val_soundspeed[iSpecies];
			k0orho		= val_normal[0] / val_density[iSpecies];
			k1orho		= val_normal[1] / val_density[iSpecies];
			if (iSpecies < nDiatomics)
				gm1 = GammaDiatomic-1.0;
			else
				gm1 = GammaMonatomic-1.0;

			gm1_o_c2	= gm1/c2;
			gm1_o_rhoxc = gm1/rhoxc;
			sqvel = val_velocity[iSpecies][0]*val_velocity[iSpecies][0]+val_velocity[iSpecies][1]*val_velocity[iSpecies][1]+val_velocity[iSpecies][2]*val_velocity[iSpecies][2];

			val_invp_tensor[loc+0][loc+0]=val_normal[0]-val_normal[2]*val_velocity[iSpecies][1] / val_density[iSpecies]+val_normal[1]*val_velocity[iSpecies][2] / val_density[iSpecies]-val_normal[0]*0.5*gm1*sqvel/c2;
			val_invp_tensor[loc+0][loc+1]=val_normal[0]*gm1*val_velocity[iSpecies][0]/c2;
			val_invp_tensor[loc+0][loc+2]=val_normal[2] / val_density[iSpecies]+val_normal[0]*gm1*val_velocity[iSpecies][1]/c2;
			val_invp_tensor[loc+0][loc+3]=-val_normal[1] / val_density[iSpecies]+val_normal[0]*gm1*val_velocity[iSpecies][2]/c2;
			val_invp_tensor[loc+0][loc+4]=-val_normal[0]*gm1/c2;

			val_invp_tensor[loc+1][loc+0]=val_normal[1]+val_normal[2]*val_velocity[iSpecies][0] / val_density[iSpecies]-val_normal[0]*val_velocity[iSpecies][2] / val_density[iSpecies]-val_normal[1]*0.5*gm1*sqvel/c2;
			val_invp_tensor[loc+1][loc+1]=-val_normal[2] / val_density[iSpecies]+val_normal[1]*gm1*val_velocity[iSpecies][0]/c2;
			val_invp_tensor[loc+1][loc+2]=val_normal[1]*gm1*val_velocity[iSpecies][1]/c2;
			val_invp_tensor[loc+1][loc+3]=val_normal[0] / val_density[iSpecies]+val_normal[1]*gm1*val_velocity[iSpecies][2]/c2;
			val_invp_tensor[loc+1][loc+4]=-val_normal[1]*gm1/c2;

			val_invp_tensor[loc+2][loc+0]=val_normal[2]-val_normal[1]*val_velocity[iSpecies][0] / val_density[iSpecies]+val_normal[0]*val_velocity[iSpecies][1] / val_density[iSpecies]-val_normal[2]*0.5*gm1*sqvel/c2;
			val_invp_tensor[loc+2][loc+1]=val_normal[1] / val_density[iSpecies]+val_normal[2]*gm1*val_velocity[iSpecies][0]/c2;
			val_invp_tensor[loc+2][loc+2]=-val_normal[0] / val_density[iSpecies]+val_normal[2]*gm1*val_velocity[iSpecies][1]/c2;
			val_invp_tensor[loc+2][loc+3]=val_normal[2]*gm1*val_velocity[iSpecies][2]/c2;
			val_invp_tensor[loc+2][loc+4]=-val_normal[2]*gm1/c2;

			val_invp_tensor[loc+3][loc+0]=-(val_normal[0]*val_velocity[iSpecies][0]+val_normal[1]*val_velocity[iSpecies][1]+val_normal[2]*val_velocity[iSpecies][2]) / val_density[iSpecies]+0.5*gm1*sqvel/rhoxc;
			val_invp_tensor[loc+3][loc+1]=val_normal[0] / val_density[iSpecies]-gm1*val_velocity[iSpecies][0]/rhoxc;
			val_invp_tensor[loc+3][loc+2]=val_normal[1] / val_density[iSpecies]-gm1*val_velocity[iSpecies][1]/rhoxc;
			val_invp_tensor[loc+3][loc+3]=val_normal[2] / val_density[iSpecies]-gm1*val_velocity[iSpecies][2]/rhoxc;
			if (iSpecies < nDiatomics)
				val_invp_tensor[loc+3][loc+4]=(GammaDiatomic-1.0)/rhoxc;
			else
				val_invp_tensor[loc+3][loc+4]=(GammaMonatomic-1.0)/rhoxc;

			val_invp_tensor[loc+4][loc+0]=(val_normal[0]*val_velocity[iSpecies][0]+val_normal[1]*val_velocity[iSpecies][1]+val_normal[2]*val_velocity[iSpecies][2]) / val_density[iSpecies]+0.5*gm1*sqvel/rhoxc;
			val_invp_tensor[loc+4][loc+1]=-val_normal[0] / val_density[iSpecies]-gm1*val_velocity[iSpecies][0]/rhoxc;
			val_invp_tensor[loc+4][loc+2]=-val_normal[1] / val_density[iSpecies]-gm1*val_velocity[iSpecies][1]/rhoxc;
			val_invp_tensor[loc+4][loc+3]=-val_normal[2] / val_density[iSpecies]-gm1*val_velocity[iSpecies][2]/rhoxc;
			if (iSpecies < nDiatomics)
				val_invp_tensor[loc+4][loc+4]=(GammaDiatomic-1.0)/rhoxc;
			else
				val_invp_tensor[loc+4][loc+4]=(GammaMonatomic-1.0)/rhoxc;
		}
	}

	if(nDim == 2) {

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

			rhoxc		= val_density[iSpecies] * val_soundspeed[iSpecies];
			c2			= val_soundspeed[iSpecies] * val_soundspeed[iSpecies];
			k0orho		= val_normal[0] / val_density[iSpecies];
			k1orho		= val_normal[1] / val_density[iSpecies];
			if (iSpecies < nDiatomics)
				gm1 = GammaDiatomic-1.0;
			else
				gm1 = GammaMonatomic-1.0;

			gm1_o_c2	= gm1/c2;
			gm1_o_rhoxc = gm1/rhoxc;

			sqvel = val_velocity[iSpecies][0]*val_velocity[iSpecies][0]+val_velocity[iSpecies][1]*val_velocity[iSpecies][1];

			val_invp_tensor[loc+0][loc+0]=1.0-0.5*gm1_o_c2*sqvel;
			val_invp_tensor[loc+0][loc+1]=gm1_o_c2*val_velocity[iSpecies][0];
			val_invp_tensor[loc+0][loc+2]=gm1_o_c2*val_velocity[iSpecies][1];
			val_invp_tensor[loc+0][loc+3]=-gm1_o_c2;

			val_invp_tensor[loc+1][loc+0]=-k1orho*val_velocity[iSpecies][0]+k0orho*val_velocity[iSpecies][1];
			val_invp_tensor[loc+1][loc+1]=k1orho;
			val_invp_tensor[loc+1][loc+2]=-k0orho;
			val_invp_tensor[loc+1][loc+3]=0.0;

			val_invp_tensor[loc+2][loc+0]=-k0orho*val_velocity[iSpecies][0]-k1orho*val_velocity[iSpecies][1]+0.5*gm1_o_rhoxc*sqvel;
			val_invp_tensor[loc+2][loc+1]=k0orho-gm1_o_rhoxc*val_velocity[iSpecies][0];
			val_invp_tensor[loc+2][loc+2]=k1orho-gm1_o_rhoxc*val_velocity[iSpecies][1];
			val_invp_tensor[loc+2][loc+3]=gm1_o_rhoxc;

			val_invp_tensor[loc+3][loc+0]=k0orho*val_velocity[iSpecies][0]+k1orho*val_velocity[iSpecies][1]+0.5*gm1_o_rhoxc*sqvel;
			val_invp_tensor[loc+3][loc+1]=-k0orho-gm1_o_rhoxc*val_velocity[iSpecies][0];
			val_invp_tensor[loc+3][loc+2]=-k1orho-gm1_o_rhoxc*val_velocity[iSpecies][1];
			val_invp_tensor[loc+3][loc+3]=gm1_o_rhoxc;
		}
	}
}

void CNumerics::GetPArtCompMatrix(double *val_density, double *val_velocity, double *val_betainc2, 
		double *val_normal, double **val_p_tensor) {
	double beta2, a, a2, velProj, area2, sx, sy, sz = 0.0, u ,v, w = 0.0;

	beta2 = *val_betainc2;
	sx = val_normal[0]; sy = val_normal[1]; 
	if (nDim == 3) sz = val_normal[2];
	u = val_velocity[0]; v = val_velocity[1]; 
	if (nDim == 3) w = val_velocity[2];
	velProj = u*sx + v*sy;
	if (nDim == 3) velProj += w*sz;
	area2 = sx*sx + sy*sy;
	if (nDim == 3) area2 += sz*sz;
	a2 = velProj*velProj + beta2*area2; a = sqrt(a2);
	double factor = 1/(2.0*beta2*a2 + EPS);

	if (nDim == 3) {
		cout << "Not implemented" << endl;
	}

	if(nDim == 2) {
		val_p_tensor[0][0] = 0.0;
		val_p_tensor[0][1] = factor*beta2*a;
		val_p_tensor[0][2] = -factor*beta2*a;

		val_p_tensor[1][0] = -factor*2.0*sy*beta2;
		val_p_tensor[1][1] = factor*(u*(a+velProj) + sx*beta2);
		val_p_tensor[1][2] = factor*(u*(velProj-a) + sx*beta2);

		val_p_tensor[2][0] = factor*2.0*sx*beta2;
		val_p_tensor[2][1] = factor*(v*(a+velProj) + sy*beta2);
		val_p_tensor[2][2] = factor*(v*(velProj-a) + sy*beta2);
	}
}

void CNumerics::GetPArtCompMatrix_inv(double *val_density, double *val_velocity, double *val_betainc2, 
		double *val_normal, double **val_invp_tensor) {
	double beta2, a, a2, velProj, area2, sx, sy, u ,v, density;

	beta2 = *val_betainc2;
	density = *val_density;

	sx = val_normal[0]; sy = val_normal[1];
	u = val_velocity[0]; v = val_velocity[1];
	velProj = u*sx + v*sy;
	area2 = sx*sx + sy*sy;
	a2 = velProj*velProj + beta2*area2; a = sqrt(a2);

	if (nDim == 3) {
		cout << "Not implemented" << endl;
	}

	if(nDim == 2) {
		val_invp_tensor[0][0] = (sy*u-sx*v)/density;
		val_invp_tensor[0][1] = -v*velProj-sy*beta2;
		val_invp_tensor[0][2] = u*velProj+sx*beta2;

		val_invp_tensor[1][0] = (a-velProj)/density;
		val_invp_tensor[1][1] = beta2*sx;
		val_invp_tensor[1][2] = beta2*sy;

		val_invp_tensor[2][0] = (-a-velProj)/density;
		val_invp_tensor[2][1] = beta2*sx;
		val_invp_tensor[2][2] = beta2*sy;

	}
}

void CNumerics::GetJacInviscidLambda_fabs(double *val_velocity, double val_soundspeed, 
		double *val_normal, double *val_Lambda_Vector) {
	double ProjVelocity = 0;

	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		ProjVelocity += val_velocity[iDim]*val_normal[iDim];

	if (nDim == 3) {
		val_Lambda_Vector[0] = fabs(ProjVelocity);
		val_Lambda_Vector[1] = fabs(ProjVelocity);
		val_Lambda_Vector[2] = fabs(ProjVelocity);
		val_Lambda_Vector[3] = fabs(ProjVelocity + val_soundspeed);
		val_Lambda_Vector[4] = fabs(ProjVelocity - val_soundspeed);
	}

	if(nDim == 2) {
		val_Lambda_Vector[0] = fabs(ProjVelocity);
		val_Lambda_Vector[1] = fabs(ProjVelocity);
		val_Lambda_Vector[2] = fabs(ProjVelocity + val_soundspeed);
		val_Lambda_Vector[3] = fabs(ProjVelocity - val_soundspeed);
	}
}

void CNumerics::ConsVar2PrimVar(double *val_consvar, double *val_primvar) {
	double Density = val_consvar[0], sq_vel = 0;
	unsigned short iDim;

	for (iDim = 0; iDim < nDim; iDim++) {
		val_primvar[iDim+1] = val_consvar[iDim+1]/Density;
		sq_vel += val_primvar[iDim+1]*val_primvar[iDim+1];
	}
	val_primvar[nVar-1] = Gamma_Minus_One*(val_consvar[nDim+1]-0.5*sq_vel*val_consvar[0]);
	val_primvar[0] = val_primvar[nVar-1] / (Gas_Constant*Density);
}

void CNumerics::ConsGrad2PrimGrad(double *val_flowsol, double **val_consvar_grad, double **val_primvar_grad) {
	unsigned short iDim, jDim;
	double grad_p, grad_u = 0.0, grad_v = 0.0, grad_w = 0.0, sqvel; 

	// Apply chain rule : Grad_Prim = (\partial(Wp)/\partial(Wc))*Grad_Cons
	sqvel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		sqvel += (val_flowsol[iDim+1]/val_flowsol[0])*(val_flowsol[iDim+1]/val_flowsol[0]);

	for (iDim = 0; iDim < nDim; iDim++){
		// Compute Gradients with respect to Primitive Variables
		grad_u = (-val_flowsol[1]/val_flowsol[0]*val_consvar_grad[0][iDim] + val_consvar_grad[1][iDim])/val_flowsol[0];  // velocity gradients.
		grad_v = (-val_flowsol[2]/val_flowsol[0]*val_consvar_grad[0][iDim] + val_consvar_grad[2][iDim])/val_flowsol[0];  
		if (nDim == 3) grad_w = (-val_flowsol[3]/val_flowsol[0]*val_consvar_grad[0][iDim] + val_consvar_grad[3][iDim])/val_flowsol[0];  

		grad_p = (Gamma-1.0)*(0.5*sqvel*val_consvar_grad[0][iDim] + val_consvar_grad[nDim+1][iDim]);   // grad_p -> pressure gradient.
		for (jDim = 0; jDim < nDim; jDim++)
			grad_p -= (Gamma-1.0)*val_flowsol[jDim+1]/val_flowsol[0]*val_consvar_grad[jDim+1][iDim];

		// Set New Gradient
		val_primvar_grad[0][iDim] = val_consvar_grad[0][iDim]; // Density
		val_primvar_grad[1][iDim] = grad_u;					   // Velocities
		val_primvar_grad[2][iDim] = grad_v;
		if (nDim == 3) val_primvar_grad[3][iDim] = grad_w;
		val_primvar_grad[nDim+1][iDim] = grad_p;               // Pressure
	}
}


void CNumerics::GetViscousFlux(double *val_primvar, double **val_gradprimvar, 
		double val_laminar_viscosity, double val_eddy_viscosity, double val_mach_inf) {

	double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	double cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
	double heat_flux_factor = cp * (val_laminar_viscosity/PRANDTL + val_eddy_viscosity/PRANDTL_TURB);

	double div_vel = 0.0;
	for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
		div_vel += val_gradprimvar[iDim+1][iDim];

	for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
		for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
			tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
																							-TWO3*total_viscosity*div_vel*delta[iDim][jDim];
		}
	}

	// Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure]
	if (nDim == 3) {
		Flux_Tensor[0][0] = 0.0; 
		Flux_Tensor[1][0] = tau[0][0]; 
		Flux_Tensor[2][0] = tau[0][1];
		Flux_Tensor[3][0] = tau[0][2];
		Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
				heat_flux_factor*val_gradprimvar[0][0];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = tau[1][0];
		Flux_Tensor[2][1] = tau[1][1]; 
		Flux_Tensor[3][1] = tau[1][2];
		Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
				heat_flux_factor*val_gradprimvar[0][1];

		Flux_Tensor[0][2] = 0.0;
		Flux_Tensor[1][2] = tau[2][0];
		Flux_Tensor[2][2] = tau[2][1];
		Flux_Tensor[3][2] = tau[2][2]; 
		Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
				heat_flux_factor*val_gradprimvar[0][2];
	}
	if (nDim == 2) {
		Flux_Tensor[0][0] = 0.0; 
		Flux_Tensor[1][0] = tau[0][0]; 
		Flux_Tensor[2][0] = tau[0][1];
		Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
				heat_flux_factor*val_gradprimvar[0][0];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = tau[1][0];
		Flux_Tensor[2][1] = tau[1][1]; 
		Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
				heat_flux_factor*val_gradprimvar[0][1];
	}
}

void CNumerics::GetViscousProjFlux(double *val_primvar, double **val_gradprimvar, double *val_normal, double val_laminar_viscosity, 
		double val_eddy_viscosity) {
	unsigned short iVar, iDim, jDim;
	double total_viscosity, heat_flux_factor, div_vel, cp; 

	total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
	heat_flux_factor = cp * (val_laminar_viscosity/PRANDTL + val_eddy_viscosity/PRANDTL_TURB);

	div_vel = 0.0;
	for (iDim = 0 ; iDim < nDim; iDim++)
		div_vel += val_gradprimvar[iDim+1][iDim];

	for (iDim = 0 ; iDim < nDim; iDim++)
		for (jDim = 0 ; jDim < nDim; jDim++)
			tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] ) - 
			TWO3*total_viscosity*div_vel*delta[iDim][jDim];

	/*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
	if (nDim == 3) {
		Flux_Tensor[0][0] = 0.0; 
		Flux_Tensor[1][0] = tau[0][0]; 
		Flux_Tensor[2][0] = tau[0][1];
		Flux_Tensor[3][0] = tau[0][2];
		Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
				heat_flux_factor*val_gradprimvar[0][0];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = tau[1][0];
		Flux_Tensor[2][1] = tau[1][1]; 
		Flux_Tensor[3][1] = tau[1][2];
		Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
				heat_flux_factor*val_gradprimvar[0][1];

		Flux_Tensor[0][2] = 0.0;
		Flux_Tensor[1][2] = tau[2][0];
		Flux_Tensor[2][2] = tau[2][1];
		Flux_Tensor[3][2] = tau[2][2]; 
		Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
				heat_flux_factor*val_gradprimvar[0][2];
	}

	if (nDim == 2) {
		Flux_Tensor[0][0] = 0.0; 
		Flux_Tensor[1][0] = tau[0][0]; 
		Flux_Tensor[2][0] = tau[0][1];
		Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
				heat_flux_factor*val_gradprimvar[0][0];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = tau[1][0];
		Flux_Tensor[2][1] = tau[1][1]; 
		Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
				heat_flux_factor*val_gradprimvar[0][1];
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Flux_Tensor[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) 
			Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
	}
}

void CNumerics::GetViscousArtCompProjFlux(double *val_primvar, double **val_gradprimvar, double *val_normal, double val_laminar_viscosity, 
		double val_eddy_viscosity) {
	unsigned short iVar, iDim, jDim;
	double total_viscosity; 

	total_viscosity = val_laminar_viscosity + val_eddy_viscosity;

	for (iDim = 0 ; iDim < nDim; iDim++)
		for (jDim = 0 ; jDim < nDim; jDim++)
			tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] );

	/*--- Gradient of primitive variables -> [Pressure vel_x vel_y vel_z] ---*/
	if (nDim == 3) {
		Flux_Tensor[0][0] = 0.0; 
		Flux_Tensor[1][0] = tau[0][0]; 
		Flux_Tensor[2][0] = tau[0][1];
		Flux_Tensor[3][0] = tau[0][2];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = tau[1][0];
		Flux_Tensor[2][1] = tau[1][1]; 
		Flux_Tensor[3][1] = tau[1][2];

		Flux_Tensor[0][2] = 0.0;
		Flux_Tensor[1][2] = tau[2][0];
		Flux_Tensor[2][2] = tau[2][1];
		Flux_Tensor[3][2] = tau[2][2]; 
	}

	if (nDim == 2) {
		Flux_Tensor[0][0] = 0.0; 
		Flux_Tensor[1][0] = tau[0][0]; 
		Flux_Tensor[2][0] = tau[0][1];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = tau[1][0];
		Flux_Tensor[2][1] = tau[1][1]; 
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Flux_Tensor[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) 
			Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
	}
}

void CNumerics::GetViscousProjJacs(double val_density, double val_pressure, double val_laminar_viscosity, 
		double val_eddy_viscosity, double val_dist_ij, double *val_normal, double val_dS,
		double *val_Mean_PrimVar, double *val_Proj_Visc_Flux, double **val_Proj_Jac_Tensor_i, double **val_Proj_Jac_Tensor_j) {

	/*--- TSL-Approximation of Viscous NS Jacobians ---*/
	unsigned short iDim, iVar, jVar;
	double theta = 0.0;
	double sqvel = 0.0;
	double proj_viscousflux_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		theta += val_normal[iDim]*val_normal[iDim];
		sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
		proj_viscousflux_vel += val_Proj_Visc_Flux[iDim+1]*val_Mean_PrimVar[iDim+1];
	}
	double phi = 0.5*(Gamma-1.)*sqvel;

	double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	double heat_flux_factor = val_laminar_viscosity / PRANDTL + val_eddy_viscosity / PRANDTL_TURB;
	double cpoR = Gamma/(Gamma-1.); // cp over R
	double factor = total_viscosity/(val_density*val_dist_ij)*val_dS;
	double phi_rho = -cpoR*heat_flux_factor*val_pressure/(val_density*val_density);
	double phi_p = cpoR*heat_flux_factor/val_density;
	double rhoovisc = val_density/total_viscosity; // rho over viscosity

	if (nDim == 3) {
		double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		double thetay = theta + val_normal[1]*val_normal[1]/3.0;
		double thetaz = theta + val_normal[2]*val_normal[2]/3.0;

		double etax = val_normal[1]*val_normal[2]/3.0;
		double etay = val_normal[0]*val_normal[2]/3.0;
		double etaz = val_normal[0]*val_normal[1]/3.0;

		double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz   + val_Mean_PrimVar[3]*etay;
		double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay + val_Mean_PrimVar[3]*etax;
		double piz = val_Mean_PrimVar[1]*etay   + val_Mean_PrimVar[2]*etax   + val_Mean_PrimVar[3]*thetaz;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[0][4] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = factor*pix;
		val_Proj_Jac_Tensor_i[1][1] = -factor*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor*etaz;
		val_Proj_Jac_Tensor_i[1][3] = -factor*etay;
		val_Proj_Jac_Tensor_i[1][4] = 0.0;
		val_Proj_Jac_Tensor_i[2][0] = factor*piy;
		val_Proj_Jac_Tensor_i[2][1] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor*thetay;
		val_Proj_Jac_Tensor_i[2][3] = -factor*etax;
		val_Proj_Jac_Tensor_i[2][4] = 0.0;
		val_Proj_Jac_Tensor_i[3][0] = factor*piz;
		val_Proj_Jac_Tensor_i[3][1] = -factor*etay;
		val_Proj_Jac_Tensor_i[3][2] = -factor*etax;
		val_Proj_Jac_Tensor_i[3][3] = -factor*thetaz;
		val_Proj_Jac_Tensor_i[3][4] = 0.0;
		val_Proj_Jac_Tensor_i[4][0] = -factor*(rhoovisc*theta*(phi_rho+phi*phi_p) - 
				pix*val_Mean_PrimVar[1]+piy*val_Mean_PrimVar[2]+piz*val_Mean_PrimVar[3]);
		val_Proj_Jac_Tensor_i[4][1] = -factor*(pix-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[1]);
		val_Proj_Jac_Tensor_i[4][2] = -factor*(piy-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[4][3] = -factor*(piz-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[3]);
		val_Proj_Jac_Tensor_i[4][4] = -factor*((Gamma-1)*rhoovisc*theta*phi_p);

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		factor = 0.5/val_density; // dS already included in Proj_Visc_Flux
		val_Proj_Jac_Tensor_i[4][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[4][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_i[4][3] += factor*val_Proj_Visc_Flux[3];
		val_Proj_Jac_Tensor_j[4][3] += factor*val_Proj_Visc_Flux[3];

	}

	if (nDim == 2) {
		double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		double thetay = theta + val_normal[1]*val_normal[1]/3.0;

		double etaz = val_normal[0]*val_normal[1]/3.0;

		double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz;
		double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = factor*pix;
		val_Proj_Jac_Tensor_i[1][1] = -factor*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor*etaz;
		val_Proj_Jac_Tensor_i[1][3] = 0.0;
		val_Proj_Jac_Tensor_i[2][0] = factor*piy;
		val_Proj_Jac_Tensor_i[2][1] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor*thetay;
		val_Proj_Jac_Tensor_i[2][3] = 0.0;
		val_Proj_Jac_Tensor_i[3][0] = -factor*(rhoovisc*theta*(phi_rho+phi*phi_p) - 
				pix*val_Mean_PrimVar[1]+piy*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[3][1] = -factor*(pix-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[1]);
		val_Proj_Jac_Tensor_i[3][2] = -factor*(piy-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[3][3] = -factor*((Gamma-1)*rhoovisc*theta*phi_p);

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		factor = 0.5/val_density; // dS already included in Proj_Visc_Flux
		val_Proj_Jac_Tensor_i[3][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[3][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[3][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[3][2] += factor*val_Proj_Visc_Flux[2];
	}
}

void CNumerics::GetViscousArtCompProjJacs(double val_density, double val_pressure, double val_laminar_viscosity, 
		double val_eddy_viscosity, double val_dist_ij, double *val_normal, double val_dS,
		double *val_Mean_PrimVar, double *val_Proj_Visc_Flux, double **val_Proj_Jac_Tensor_i, double **val_Proj_Jac_Tensor_j) {

	/*--- TSL-Approximation of Viscous NS Jacobians ---*/
	unsigned short iDim, iVar, jVar;
	double theta = 0.0;
	double sqvel = 0.0;
	double proj_viscousflux_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		theta += val_normal[iDim]*val_normal[iDim];
		sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
		proj_viscousflux_vel += val_Proj_Visc_Flux[iDim+1]*val_Mean_PrimVar[iDim+1];
	}
	double phi = 0.5*(Gamma-1.)*sqvel;

	double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	double heat_flux_factor = val_laminar_viscosity / PRANDTL + val_eddy_viscosity / PRANDTL_TURB;
	double cpoR = Gamma/(Gamma-1.); // cp over R
	double factor = total_viscosity/(val_density*val_dist_ij)*val_dS;
	double phi_rho = -cpoR*heat_flux_factor*val_pressure/(val_density*val_density);
	double phi_p = cpoR*heat_flux_factor/val_density;
	double rhoovisc = val_density/total_viscosity; // rho over viscosity

	if (nDim == 3) {
		double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		double thetay = theta + val_normal[1]*val_normal[1]/3.0;
		double thetaz = theta + val_normal[2]*val_normal[2]/3.0;

		double etax = val_normal[1]*val_normal[2]/3.0;
		double etay = val_normal[0]*val_normal[2]/3.0;
		double etaz = val_normal[0]*val_normal[1]/3.0;

		double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz   + val_Mean_PrimVar[3]*etay;
		double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay + val_Mean_PrimVar[3]*etax;
		double piz = val_Mean_PrimVar[1]*etay   + val_Mean_PrimVar[2]*etax   + val_Mean_PrimVar[3]*thetaz;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[0][4] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = factor*pix;
		val_Proj_Jac_Tensor_i[1][1] = -factor*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor*etaz;
		val_Proj_Jac_Tensor_i[1][3] = -factor*etay;
		val_Proj_Jac_Tensor_i[1][4] = 0.0;
		val_Proj_Jac_Tensor_i[2][0] = factor*piy;
		val_Proj_Jac_Tensor_i[2][1] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor*thetay;
		val_Proj_Jac_Tensor_i[2][3] = -factor*etax;
		val_Proj_Jac_Tensor_i[2][4] = 0.0;
		val_Proj_Jac_Tensor_i[3][0] = factor*piz;
		val_Proj_Jac_Tensor_i[3][1] = -factor*etay;
		val_Proj_Jac_Tensor_i[3][2] = -factor*etax;
		val_Proj_Jac_Tensor_i[3][3] = -factor*thetaz;
		val_Proj_Jac_Tensor_i[3][4] = 0.0;
		val_Proj_Jac_Tensor_i[4][0] = -factor*(rhoovisc*theta*(phi_rho+phi*phi_p) - 
				pix*val_Mean_PrimVar[1]+piy*val_Mean_PrimVar[2]+piz*val_Mean_PrimVar[3]);
		val_Proj_Jac_Tensor_i[4][1] = -factor*(pix-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[1]);
		val_Proj_Jac_Tensor_i[4][2] = -factor*(piy-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[4][3] = -factor*(piz-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[3]);
		val_Proj_Jac_Tensor_i[4][4] = -factor*((Gamma-1)*rhoovisc*theta*phi_p);

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		factor = 0.5/val_density; // dS already included in Proj_Visc_Flux
		val_Proj_Jac_Tensor_i[4][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[4][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_i[4][3] += factor*val_Proj_Visc_Flux[3];
		val_Proj_Jac_Tensor_j[4][3] += factor*val_Proj_Visc_Flux[3];

	}

	if (nDim == 2) {
		double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		double thetay = theta + val_normal[1]*val_normal[1]/3.0;

		double etaz = val_normal[0]*val_normal[1]/3.0;

		double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz;
		double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = factor*pix;
		val_Proj_Jac_Tensor_i[1][1] = -factor*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor*etaz;
		val_Proj_Jac_Tensor_i[1][3] = 0.0;
		val_Proj_Jac_Tensor_i[2][0] = factor*piy;
		val_Proj_Jac_Tensor_i[2][1] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor*thetay;
		val_Proj_Jac_Tensor_i[2][3] = 0.0;
		val_Proj_Jac_Tensor_i[3][0] = -factor*(rhoovisc*theta*(phi_rho+phi*phi_p) - 
				pix*val_Mean_PrimVar[1]+piy*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[3][1] = -factor*(pix-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[1]);
		val_Proj_Jac_Tensor_i[3][2] = -factor*(piy-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[3][3] = -factor*((Gamma-1)*rhoovisc*theta*phi_p);

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		factor = 0.5/val_density; // dS already included in Proj_Visc_Flux
		val_Proj_Jac_Tensor_i[3][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[3][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[3][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[3][2] += factor*val_Proj_Visc_Flux[2];
	}
}
