/*!
 * \file numerics_structure.cpp
 * \brief This file contains all the numerical methods.
 * \author F. Palacios, T. Economon
 * \version 4.2.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/numerics_structure.hpp"

CNumerics::CNumerics(void) {

  Normal = NULL;
  UnitNormal = NULL;
  U_n = NULL;
  U_nM1 = NULL;
  U_nP1  = NULL;
  Proj_Flux_Tensor  = NULL;
  Flux_Tensor = NULL;
  tau  = NULL;
  delta  = NULL;

  Diffusion_Coeff_i = NULL;
  Diffusion_Coeff_j = NULL;

  Enthalpy_formation = NULL;
  Theta_v = NULL;
  var = NULL;

  Ys          = NULL;
  dFdYj       = NULL;
  dFdYi       = NULL;
  sumdFdYih   = NULL;
  sumdFdYjh   = NULL;
  sumdFdYieve = NULL;
  sumdFdYjeve = NULL;
}

CNumerics::CNumerics(unsigned short val_nDim, unsigned short val_nVar,
                     CConfig *config) {
  
  unsigned short iVar, iDim, jDim;
  
	nDim = val_nDim;
	nVar = val_nVar;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  Prandtl_Lam = config->GetPrandtl_Lam();
  Prandtl_Turb = config->GetPrandtl_Turb();
	Gas_Constant = config->GetGas_ConstantND();

	UnitNormal = new su2double [nDim];
	UnitNormald = new su2double [nDim];

	Normal = new su2double [nDim];
	Flux_Tensor = new su2double* [nVar];
	for (iVar = 0; iVar < (nVar); iVar++)
		Flux_Tensor[iVar] = new su2double [nDim];

	tau = new su2double* [nDim];
	delta = new su2double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		tau[iDim] = new su2double [nDim];
		delta[iDim] = new su2double [nDim];
	}

	for (iDim = 0; iDim < nDim; iDim++) {
		for (jDim = 0; jDim < nDim; jDim++) {
			if (iDim == jDim) delta[iDim][jDim]=1.0;
			else delta[iDim][jDim]=0.0;
		}
	}

	U_n = new su2double [nVar];
	U_nM1 = new su2double [nVar];
	U_nP1 = new su2double [nVar];

	Proj_Flux_Tensor = new su2double [nVar];

	turb_ke_i = 0.0;
	turb_ke_j = 0.0;
  
  Diffusion_Coeff_i = NULL;
  Diffusion_Coeff_j = NULL;
  
  Enthalpy_formation = NULL;
  Theta_v = NULL;
  var = NULL;

  Ys          = NULL;
  dFdYj       = NULL;
  dFdYi       = NULL;
  sumdFdYih   = NULL;
  sumdFdYjh   = NULL;
  sumdFdYieve = NULL;
  sumdFdYjeve = NULL;
    
  Vector = new su2double[nDim];
  
  l = new su2double [nDim];
  m = new su2double [nDim];
  
}

CNumerics::~CNumerics(void) {
  if (Normal!=NULL) delete [] Normal;
  if (UnitNormal!=NULL) delete [] UnitNormal;

  if (U_n!=NULL) delete [] U_n;
  if (U_nM1!=NULL) delete [] U_nM1;
  if (U_nP1!=NULL) delete [] U_nP1;

	// visc
  if (Proj_Flux_Tensor!=NULL) delete [] Proj_Flux_Tensor;

  if (Flux_Tensor!=NULL){
    for (unsigned short iVar = 0; iVar < nDim+3; iVar++) {
      delete [] Flux_Tensor[iVar];
    }
    delete [] Flux_Tensor;
  }

  if (tau!=NULL && delta != NULL){
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] tau[iDim];
      delete [] delta[iDim];
    }
    delete [] tau;
    delete [] delta;
  }

	if (Ys != NULL) delete [] Ys;
  if (sumdFdYih != NULL) delete [] sumdFdYih;
  if (sumdFdYjh != NULL) delete [] sumdFdYjh;
  if (sumdFdYieve != NULL) delete [] sumdFdYieve;
  if (sumdFdYjeve != NULL) delete [] sumdFdYjeve;
  if (Diffusion_Coeff_i != NULL) delete [] Diffusion_Coeff_i;
  if (Diffusion_Coeff_j != NULL) delete [] Diffusion_Coeff_j;
  if (Vector != NULL) delete [] Vector;
  if (var != NULL) delete [] var;

	if(Enthalpy_formation != NULL) delete [] Enthalpy_formation;
	if(Theta_v != NULL) delete [] Theta_v;

}

void CNumerics::GetInviscidFlux(su2double val_density, su2double *val_velocity,
		su2double val_pressure, su2double val_enthalpy) {
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
	if (nDim == 2) {
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

void CNumerics::GetInviscidProjFlux(su2double *val_density,
                                    su2double *val_velocity,
                                    su2double *val_pressure,
                                    su2double *val_enthalpy,
                                    su2double *val_normal,
                                    su2double *val_Proj_Flux) {
  
    su2double rhou, rhov, rhow;
    
	if (nDim == 2) {
    
		rhou = (*val_density)*val_velocity[0];
		rhov = (*val_density)*val_velocity[1];

		val_Proj_Flux[0] = rhou*val_normal[0];
		val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
		val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
		val_Proj_Flux[3] = rhou*(*val_enthalpy)*val_normal[0];

		val_Proj_Flux[0] += rhov*val_normal[1];
		val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
		val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
		val_Proj_Flux[3] += rhov*(*val_enthalpy)*val_normal[1];
    
	} 
	else {
    
		rhou = (*val_density)*val_velocity[0];
		rhov = (*val_density)*val_velocity[1];
		rhow = (*val_density)*val_velocity[2];

		val_Proj_Flux[0] = rhou*val_normal[0];
		val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
		val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
		val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
		val_Proj_Flux[4] = rhou*(*val_enthalpy)*val_normal[0];

		val_Proj_Flux[0] += rhov*val_normal[1];
		val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
		val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
		val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
		val_Proj_Flux[4] += rhov*(*val_enthalpy)*val_normal[1];

		val_Proj_Flux[0] += rhow*val_normal[2];
		val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
		val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
		val_Proj_Flux[3] += (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
		val_Proj_Flux[4] += rhow*(*val_enthalpy)*val_normal[2];
    
	}

}

void CNumerics::GetInviscidArtCompProjFlux(su2double *val_density,
                                           su2double *val_velocity,
                                           su2double *val_pressure,
                                           su2double *val_betainc2,
                                           su2double *val_normal,
                                           su2double *val_Proj_Flux) {
    su2double rhou, rhov, rhow;
    
   	if (nDim == 2) {
      rhou = (*val_density)*val_velocity[0];
      rhov = (*val_density)*val_velocity[1];
      
      val_Proj_Flux[0] = (*val_betainc2)*(val_velocity[0]*val_normal[0] + val_velocity[1]*val_normal[1]);
      val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0] + rhou*val_velocity[1]*val_normal[1];
      val_Proj_Flux[2] = rhov*val_velocity[0]*val_normal[0] + (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
	}
  else {
    rhou = (*val_density)*val_velocity[0];
		rhov = (*val_density)*val_velocity[1];
		rhow = (*val_density)*val_velocity[2];
    
		val_Proj_Flux[0] = (*val_betainc2)*(val_velocity[0]*val_normal[0] + val_velocity[1]*val_normal[1] + val_velocity[2]*val_normal[2]);
		val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0] + rhou*val_velocity[1]*val_normal[1] + rhou*val_velocity[2]*val_normal[2];
		val_Proj_Flux[2] = rhov*val_velocity[0]*val_normal[0] + (rhov*val_velocity[1]+(*val_pressure))*val_normal[1] + rhov*val_velocity[2]*val_normal[2];
		val_Proj_Flux[3] = rhow*val_velocity[0]*val_normal[0] + rhow*val_velocity[1]*val_normal[1] + (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
	}
  
}

void CNumerics::GetInviscidArtComp_FreeSurf_ProjFlux(su2double *val_density, su2double *val_velocity, su2double *val_pressure, su2double *val_betainc2,
                                                     su2double *val_levelset, su2double *val_normal, su2double *val_Proj_Flux) {
  
  su2double rhou, rhov, rhow, ProjVel;
  
  if (nDim == 2) {
    rhou = (*val_density)*val_velocity[0];
		rhov = (*val_density)*val_velocity[1];
    ProjVel = (val_velocity[0]*val_normal[0] + val_velocity[1]*val_normal[1]);
    
		val_Proj_Flux[0] = (*val_betainc2)*ProjVel;
		val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0] + rhou*val_velocity[1]*val_normal[1];
		val_Proj_Flux[2] = rhov*val_velocity[0]*val_normal[0] + (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] = (*val_levelset)*ProjVel;
	}
  else {
    rhou = (*val_density)*val_velocity[0];
		rhov = (*val_density)*val_velocity[1];
		rhow = (*val_density)*val_velocity[2];
    ProjVel = (val_velocity[0]*val_normal[0] + val_velocity[1]*val_normal[1] + val_velocity[2]*val_normal[2]);

		val_Proj_Flux[0] = (*val_betainc2)*ProjVel;
		val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0] + rhou*val_velocity[1]*val_normal[1] + rhou*val_velocity[2]*val_normal[2];
		val_Proj_Flux[2] = rhov*val_velocity[0]*val_normal[0] + (rhov*val_velocity[1]+(*val_pressure))*val_normal[1] + rhov*val_velocity[2]*val_normal[2];
		val_Proj_Flux[3] = rhow*val_velocity[0]*val_normal[0] + rhow*val_velocity[1]*val_normal[1] + (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] = (*val_levelset)*ProjVel;
	}
}

void CNumerics::GetInviscidProjJac(su2double *val_velocity, su2double *val_energy,
                                   su2double *val_normal, su2double val_scale,
                                   su2double **val_Proj_Jac_Tensor) {
  AD_BEGIN_PASSIVE
  unsigned short iDim, jDim;
  su2double sqvel, proj_vel, phi, a1, a2;
  
  sqvel = 0.0, proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel    += val_velocity[iDim]*val_velocity[iDim];
    proj_vel += val_velocity[iDim]*val_normal[iDim];
  }
  
  phi = 0.5*Gamma_Minus_One*sqvel;
  a1 = Gamma*(*val_energy)-phi;
  a2 = Gamma-1.0;
  
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
  AD_END_PASSIVE
}


void CNumerics::GetInviscidProjJac(su2double *val_velocity, su2double *val_enthalpy,
		su2double *val_chi, su2double *val_kappa,
		su2double *val_normal, su2double val_scale,
		su2double **val_Proj_Jac_Tensor) {
  AD_BEGIN_PASSIVE
	unsigned short iDim, jDim;
	su2double sqvel, proj_vel, phi, a1, a2;

	sqvel = 0.0, proj_vel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		sqvel += val_velocity[iDim]*val_velocity[iDim];
		proj_vel += val_velocity[iDim]*val_normal[iDim];
	}

	phi = *val_chi + 0.5*sqvel*(*val_kappa);
	a1 = *val_enthalpy;
	a2 = *val_kappa;

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
	val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*(a2+1)*proj_vel;
  AD_END_PASSIVE
}

void CNumerics::GetInviscidArtCompProjJac(su2double *val_density, su2double *val_velocity, su2double *val_betainc2, su2double *val_normal,
		su2double val_scale, su2double **val_Proj_Jac_Tensor) {
  AD_BEGIN_PASSIVE
	unsigned short iDim;
	su2double proj_vel;

	proj_vel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		proj_vel += val_velocity[iDim]*val_normal[iDim];

    if (nDim == 2) {
		val_Proj_Jac_Tensor[0][0] = 0.0;
		val_Proj_Jac_Tensor[0][1] = val_scale*(*val_betainc2)*val_normal[0]/(*val_density);
		val_Proj_Jac_Tensor[0][2] = val_scale*(*val_betainc2)*val_normal[1]/(*val_density);
        
		val_Proj_Jac_Tensor[1][0] = val_scale*val_normal[0];
		val_Proj_Jac_Tensor[1][1] = val_scale*(val_velocity[0]*val_normal[0] + proj_vel);
		val_Proj_Jac_Tensor[1][2] = val_scale*val_velocity[0]*val_normal[1];
        
		val_Proj_Jac_Tensor[2][0] = val_scale*val_normal[1];
		val_Proj_Jac_Tensor[2][1] = val_scale*val_velocity[1]*val_normal[0];
		val_Proj_Jac_Tensor[2][2] = val_scale*(val_velocity[1]*val_normal[1] + proj_vel);
	}
	else {
		val_Proj_Jac_Tensor[0][0] = 0.0;
		val_Proj_Jac_Tensor[0][1] = val_scale*(*val_betainc2)*val_normal[0]/(*val_density);
		val_Proj_Jac_Tensor[0][2] = val_scale*(*val_betainc2)*val_normal[1]/(*val_density);
		val_Proj_Jac_Tensor[0][3] = val_scale*(*val_betainc2)*val_normal[2]/(*val_density);

		val_Proj_Jac_Tensor[1][0] = val_scale*val_normal[0];
		val_Proj_Jac_Tensor[1][1] = val_scale*(val_velocity[0]*val_normal[0] + proj_vel);
		val_Proj_Jac_Tensor[1][2] = val_scale*val_velocity[0]*val_normal[1];
		val_Proj_Jac_Tensor[1][3] = val_scale*val_velocity[0]*val_normal[2];

		val_Proj_Jac_Tensor[2][0] = val_scale*val_normal[1];
		val_Proj_Jac_Tensor[2][1] = val_scale*val_velocity[1]*val_normal[0];
		val_Proj_Jac_Tensor[2][2] = val_scale*(val_velocity[1]*val_normal[1] + proj_vel);
		val_Proj_Jac_Tensor[2][3] = val_scale*val_velocity[1]*val_normal[2];

		val_Proj_Jac_Tensor[3][0] = val_scale*val_normal[2];
		val_Proj_Jac_Tensor[3][1] = val_scale*val_velocity[2]*val_normal[0];
		val_Proj_Jac_Tensor[3][2] = val_scale*val_velocity[2]*val_normal[1];
		val_Proj_Jac_Tensor[3][3] = val_scale*(val_velocity[2]*val_normal[2] + proj_vel);
	}
  AD_END_PASSIVE
}

void CNumerics::GetInviscidArtComp_FreeSurf_ProjJac(su2double *val_density, su2double *val_ddensity, su2double *val_velocity, su2double *val_betainc2, su2double *val_levelset, su2double *val_normal,
                                                    su2double val_scale, su2double **val_Proj_Jac_Tensor) {
  
	su2double a = 0.0, b = 0.0, c = 0.0, d = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, u = 0.0, v = 0.0, w = 0.0;
  
  a = (*val_betainc2)/(*val_density);
  b = (*val_levelset)/(*val_density);
  c = (*val_ddensity);
  
  if (nDim == 2) {
    
    nx = val_normal[0];   ny = val_normal[1];
    u = val_velocity[0];  v = val_velocity[1];  d = u*nx + v*ny;
    
    val_Proj_Jac_Tensor[0][0] = 0.0;
    val_Proj_Jac_Tensor[0][1] = val_scale*a*nx;
    val_Proj_Jac_Tensor[0][2] = val_scale*a*ny;
    val_Proj_Jac_Tensor[0][3] = - val_scale*a*c*d;
    
    
    val_Proj_Jac_Tensor[1][0] = val_scale*nx;
    val_Proj_Jac_Tensor[1][1] = val_scale*(d + nx*u);
    val_Proj_Jac_Tensor[1][2] = val_scale*ny*u;
    val_Proj_Jac_Tensor[1][3] = -val_scale*c*d*u;
    
    
    val_Proj_Jac_Tensor[2][0] = val_scale*ny;
    val_Proj_Jac_Tensor[2][1] = val_scale*nx*v;
    val_Proj_Jac_Tensor[2][2] = val_scale*(d + ny*v);
    val_Proj_Jac_Tensor[2][3] = -val_scale*c*d*v;
    
    
    val_Proj_Jac_Tensor[3][0] = 0.0;
    val_Proj_Jac_Tensor[3][1] = val_scale*b*nx;
    val_Proj_Jac_Tensor[3][2] = val_scale*b*ny;
    val_Proj_Jac_Tensor[3][3] = val_scale*(d - b*c*d);
    
  }
	else {
  
    nx = val_normal[0];   ny = val_normal[1];   nz = val_normal[2];
    u = val_velocity[0];  v = val_velocity[1];  w = val_velocity[2];  d = u*nx + v*ny + w*nz;
  
    val_Proj_Jac_Tensor[0][0] = 0.0;
    val_Proj_Jac_Tensor[0][1] = val_scale*a*nx;
    val_Proj_Jac_Tensor[0][2] = val_scale*a*ny;
    val_Proj_Jac_Tensor[0][3] = val_scale*a*nz;
    val_Proj_Jac_Tensor[0][4] = - val_scale*a*c*d;
    
    
    val_Proj_Jac_Tensor[1][0] = val_scale*nx;
    val_Proj_Jac_Tensor[1][1] = val_scale*(d + nx*u);
    val_Proj_Jac_Tensor[1][2] = val_scale*ny*u;
    val_Proj_Jac_Tensor[1][3] = val_scale*nz*u;
    val_Proj_Jac_Tensor[1][4] = -val_scale*c*d*u;
    
    
    val_Proj_Jac_Tensor[2][0] = val_scale*ny;
    val_Proj_Jac_Tensor[2][1] = val_scale*nx*v;
    val_Proj_Jac_Tensor[2][2] = val_scale*(d + ny*v);
    val_Proj_Jac_Tensor[2][3] = val_scale*nz*v;
    val_Proj_Jac_Tensor[2][4] = -val_scale*c*d*v;
    
    val_Proj_Jac_Tensor[3][0] = val_scale*nz;
    val_Proj_Jac_Tensor[3][1] = val_scale*nx*w;
    val_Proj_Jac_Tensor[3][2] = val_scale*ny*w;
    val_Proj_Jac_Tensor[3][3] = val_scale*(d + nz*w);
    val_Proj_Jac_Tensor[3][4] = -val_scale*c*d*w;
    
    val_Proj_Jac_Tensor[3][0] = 0.0;
    val_Proj_Jac_Tensor[3][1] = val_scale*b*nx;
    val_Proj_Jac_Tensor[3][2] = val_scale*b*ny;
    val_Proj_Jac_Tensor[3][3] = val_scale*b*nz;
    val_Proj_Jac_Tensor[3][4] = val_scale*(d - b*c*d);
    
  }
  
}

void CNumerics::SetPastSol (su2double *val_u_nM1, su2double *val_u_n, su2double *val_u_nP1) {
	unsigned short iVar;

	for (iVar = 0; iVar < nVar; iVar++) {
		U_nM1[iVar] = val_u_nM1[iVar];
		U_n[iVar] = val_u_n[iVar];
		U_nP1[iVar] = val_u_nP1[iVar];
	}

}
void CNumerics::SetPastVolume (su2double val_volume_nM1, su2double val_volume_n, su2double val_volume_nP1) {
	Volume_nM1 = val_volume_nM1;
	Volume_n = val_volume_n;
	Volume_nP1 = val_volume_nP1;
}


void CNumerics::GetPMatrix(su2double *val_density, su2double *val_velocity,
                           su2double *val_soundspeed, su2double *val_normal, su2double **val_p_tensor) {
  
  su2double sqvel, rhooc, rhoxc;
  //su2double c2;
  
  rhooc = *val_density / *val_soundspeed;
  rhoxc = *val_density * *val_soundspeed;
  //c2 = *val_soundspeed * *val_soundspeed;
  
  if (nDim == 2) {
    
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
  else {
    
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
  
}

void CNumerics::GetPMatrix(su2double *val_density, su2double *val_velocity,
		su2double *val_soundspeed, su2double *val_enthalpy, su2double *val_chi, su2double *val_kappa, su2double *val_normal, su2double **val_p_tensor) {

  su2double sqvel, rhooc, zeta;
  //su2double rhoxc, c2;
  
  rhooc = *val_density / *val_soundspeed;
  //rhoxc = *val_density * *val_soundspeed;
  //c2 = *val_soundspeed * *val_soundspeed;

	if (nDim == 2) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];
		zeta = sqvel - (*val_kappa*0.5*sqvel + *val_chi)/(*val_kappa);

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

		val_p_tensor[3][0]= zeta;
		val_p_tensor[3][1]=*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
		val_p_tensor[3][2]=0.5*(*val_enthalpy*rhooc+*val_density*val_velocity[0]*val_normal[0]+*val_density*val_velocity[1]*val_normal[1]);
		val_p_tensor[3][3]=0.5*(*val_enthalpy*rhooc-*val_density*val_velocity[0]*val_normal[0]-*val_density*val_velocity[1]*val_normal[1]);
	}
	else {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
		zeta = sqvel - (*val_kappa*0.5*sqvel + *val_chi)/(*val_kappa);

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

		val_p_tensor[4][0]=zeta*val_normal[0]+*val_density*val_velocity[1]*val_normal[2]-*val_density*val_velocity[2]*val_normal[1];
		val_p_tensor[4][1]=zeta*val_normal[1]-*val_density*val_velocity[0]*val_normal[2]+*val_density*val_velocity[2]*val_normal[0];
		val_p_tensor[4][2]=zeta*val_normal[2]+*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
		val_p_tensor[4][3]=0.5*(*val_enthalpy*rhooc+*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2]));
		val_p_tensor[4][4]=0.5*(*val_enthalpy*rhooc-*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2]));
	}

}

void CNumerics::GetPMatrix_inv(su2double *val_density, su2double *val_velocity,
		su2double *val_soundspeed, su2double *val_normal, su2double **val_invp_tensor) {
  
	su2double rhoxc, c2, gm1, k0orho, k1orho, gm1_o_c2, gm1_o_rhoxc, sqvel;

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
	if (nDim == 2) {
    
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

void CNumerics::GetPMatrix_inv(su2double **val_invp_tensor, su2double *val_density, su2double *val_velocity,
		su2double *val_soundspeed, su2double *val_chi, su2double *val_kappa, su2double *val_normal) {

	su2double rhoxc, c2, k0orho, k1orho, sqvel, k_o_c2, k_o_rhoxc, dp_drho;

	rhoxc = *val_density * *val_soundspeed;
	c2 = *val_soundspeed * *val_soundspeed;
	k0orho = val_normal[0] / *val_density;
	k1orho = val_normal[1] / *val_density;
	k_o_c2 = (*val_kappa)/c2;
	k_o_rhoxc = (*val_kappa)/rhoxc;


	if (nDim == 3) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
		dp_drho = *val_chi + 0.5*sqvel*(*val_kappa);

		val_invp_tensor[0][0]=val_normal[0]-val_normal[2]*val_velocity[1] / *val_density + val_normal[1]*val_velocity[2] / *val_density - val_normal[0]*dp_drho/c2;
		val_invp_tensor[0][1]=val_normal[0]*val_velocity[0]*k_o_c2;
		val_invp_tensor[0][2]=val_normal[2] / *val_density + val_normal[0]*val_velocity[1]*k_o_c2;
		val_invp_tensor[0][3]=-val_normal[1] / *val_density + val_normal[0]*val_velocity[2]*k_o_c2;
		val_invp_tensor[0][4]=-val_normal[0]*k_o_c2;

		val_invp_tensor[1][0]=val_normal[1]+val_normal[2]*val_velocity[0] / *val_density - val_normal[0]*val_velocity[2] / *val_density - val_normal[1]*dp_drho/c2;
		val_invp_tensor[1][1]=-val_normal[2] / *val_density + val_normal[1]*val_velocity[0]*k_o_c2;
		val_invp_tensor[1][2]=val_normal[1]*val_velocity[1]*k_o_c2;
		val_invp_tensor[1][3]=val_normal[0] / *val_density + val_normal[1]*val_velocity[2]*k_o_c2;
		val_invp_tensor[1][4]=-val_normal[1]*k_o_c2;

		val_invp_tensor[2][0]=val_normal[2]-val_normal[1]*val_velocity[0] / *val_density + val_normal[0]*val_velocity[1] / *val_density - val_normal[2]*dp_drho/c2;
		val_invp_tensor[2][1]=val_normal[1] / *val_density + val_normal[2]*val_velocity[0]*k_o_c2;
		val_invp_tensor[2][2]=-val_normal[0] / *val_density + val_normal[2]*val_velocity[1]*k_o_c2;
		val_invp_tensor[2][3]=val_normal[2]*val_velocity[2]*k_o_c2;
		val_invp_tensor[2][4]=-val_normal[2]*k_o_c2;

		val_invp_tensor[3][0]=-(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+ dp_drho/rhoxc;
		val_invp_tensor[3][1]=val_normal[0] / *val_density - val_velocity[0]*k_o_rhoxc;
		val_invp_tensor[3][2]=val_normal[1] / *val_density- val_velocity[1]*k_o_rhoxc;
		val_invp_tensor[3][3]=val_normal[2] / *val_density- val_velocity[2]*k_o_rhoxc;
		val_invp_tensor[3][4]= k_o_rhoxc;

		val_invp_tensor[4][0]=(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+ dp_drho/rhoxc;
		val_invp_tensor[4][1]=-val_normal[0] / *val_density- val_velocity[0]*k_o_rhoxc;
		val_invp_tensor[4][2]=-val_normal[1] / *val_density- val_velocity[1]*k_o_rhoxc;
		val_invp_tensor[4][3]=-val_normal[2] / *val_density- val_velocity[2]*k_o_rhoxc;
		val_invp_tensor[4][4]= k_o_rhoxc;
	}
	if (nDim == 2) {
		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];
		dp_drho = *val_chi + 0.5*sqvel*(*val_kappa);

		val_invp_tensor[0][0]=1.0 - dp_drho/c2;
		val_invp_tensor[0][1]= k_o_c2*val_velocity[0];
		val_invp_tensor[0][2]= k_o_c2*val_velocity[1];
		val_invp_tensor[0][3]=-k_o_c2;

		val_invp_tensor[1][0]=-k1orho*val_velocity[0]+k0orho*val_velocity[1];
		val_invp_tensor[1][1]=k1orho;
		val_invp_tensor[1][2]=-k0orho;
		val_invp_tensor[1][3]=0.0;

		val_invp_tensor[2][0]=-k0orho*val_velocity[0]-k1orho*val_velocity[1] + dp_drho/rhoxc;
		val_invp_tensor[2][1]=k0orho - k_o_rhoxc*val_velocity[0];
		val_invp_tensor[2][2]=k1orho - k_o_rhoxc*val_velocity[1];
		val_invp_tensor[2][3]=k_o_rhoxc;

		val_invp_tensor[3][0]=k0orho*val_velocity[0]+k1orho*val_velocity[1] + dp_drho/rhoxc;
		val_invp_tensor[3][1]=-k0orho - k_o_rhoxc*val_velocity[0];
		val_invp_tensor[3][2]=-k1orho - k_o_rhoxc*val_velocity[1];
		val_invp_tensor[3][3]= k_o_rhoxc;
	}
}

void CNumerics::GetinvRinvPe(su2double Beta2, su2double val_enthalpy,
                             su2double val_soundspeed, su2double val_density,
                             su2double* val_velocity, su2double **invRinvPe) {

	su2double sqvel;
	su2double factor = 1.0/(val_soundspeed*val_soundspeed*Beta2);

	if (nDim == 2) {

		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

		invRinvPe[0][0] = factor;
		invRinvPe[0][1] = 0.0;
		invRinvPe[0][2] = 0.0;
		invRinvPe[0][3] = -val_density/Gamma;

		invRinvPe[1][0] = val_velocity[0]*factor;
		invRinvPe[1][1] = val_density;
		invRinvPe[1][2] = 0.0;
		invRinvPe[1][3] = -val_density*val_velocity[0]/Gamma;

		invRinvPe[2][0] = val_velocity[1]*factor;
		invRinvPe[2][1] = 0;
		invRinvPe[2][2] = val_density;
		invRinvPe[2][3] = -val_density*val_velocity[1]/Gamma;

		invRinvPe[3][0] = val_enthalpy*factor;
		invRinvPe[3][1] = val_density*val_velocity[0];
		invRinvPe[3][2] = val_density*val_velocity[1];
		invRinvPe[3][3] = -val_density*sqvel/(2.0*Gamma);
	}
	else {

		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

		invRinvPe[0][0] =  factor;
		invRinvPe[0][1] = 0.0;
		invRinvPe[0][2] = 0.0;
		invRinvPe[0][3] = 0.0;
		invRinvPe[0][4] = -val_density/Gamma;

		invRinvPe[1][0] = val_velocity[0]*factor;
		invRinvPe[1][1] = val_density;
		invRinvPe[1][2] = 0.0;
		invRinvPe[1][3] = 0.0;
		invRinvPe[1][4] = -val_density*val_velocity[0]/Gamma;

		invRinvPe[2][0] = val_velocity[1]*factor;
		invRinvPe[2][1] = 0;
		invRinvPe[2][2] = val_density;
		invRinvPe[2][3] = 0.0;
		invRinvPe[2][4] = -val_density*val_velocity[1]/Gamma;


		invRinvPe[3][0] = val_velocity[2]*factor;
		invRinvPe[3][1] = 0;
		invRinvPe[3][2] = 0;
		invRinvPe[3][3] = val_density;
		invRinvPe[3][4] = -val_density*val_velocity[2]/Gamma;

		invRinvPe[4][0] = val_enthalpy*factor;
		invRinvPe[4][1] = val_density*val_velocity[0];
		invRinvPe[4][2] = val_density*val_velocity[1];
		invRinvPe[4][3] = val_density*val_velocity[2];
		invRinvPe[4][4] = -val_density*sqvel/(2.0*Gamma);

	}

}

void CNumerics::GetRMatrix(su2double val_pressure, su2double val_soundspeed, su2double val_density, su2double* val_velocity, su2double **R_Matrix) {

	su2double sqvel;
	//su2double factor = 1.0/(val_soundspeed*val_soundspeed*1);
	su2double gm1 = Gamma - 1.0;

	if (nDim == 2) {

		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

		R_Matrix[0][0] =  0.5*gm1*sqvel;
		R_Matrix[0][1] = -val_velocity[0]*gm1;
		R_Matrix[0][2] = -val_velocity[1]*gm1;
		R_Matrix[0][3] = gm1;

		R_Matrix[1][0] = -val_velocity[0]/val_density;
		R_Matrix[1][1] = 1.0/val_density;
		R_Matrix[1][2] = 0.0;
		R_Matrix[1][3] = 0.0;

		R_Matrix[2][0] = -val_velocity[1]/val_density;
		R_Matrix[2][1] = 0.0;
		R_Matrix[2][2] = 1.0/val_density;
		R_Matrix[2][3] = 0.0;

		R_Matrix[3][0] = 0.5*gm1*sqvel/val_pressure - Gamma/val_density;
		R_Matrix[3][1] = -gm1*val_velocity[0]/val_pressure;
		R_Matrix[3][2] = -gm1*val_velocity[1]/val_pressure;
		R_Matrix[3][3] = gm1/val_pressure;
	}
	else {

		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

		R_Matrix[0][0] =  0.5*gm1*sqvel;
		R_Matrix[0][1] = -val_velocity[0]*gm1;
		R_Matrix[0][2] = -val_velocity[1]*gm1;
		R_Matrix[0][3] = -val_velocity[2]*gm1;
		R_Matrix[0][4] = gm1;

		R_Matrix[1][0] = -val_velocity[0]/val_density;
		R_Matrix[1][1] = 1.0/val_density;
		R_Matrix[1][2] = 0.0;
		R_Matrix[1][3] = 0.0;
		R_Matrix[1][4] = 0.0;

		R_Matrix[2][0] = -val_velocity[1]/val_density;
		R_Matrix[2][1] = 0.0;
		R_Matrix[2][2] = 1.0/val_density;
		R_Matrix[2][3] = 0.0;
		R_Matrix[2][4] = 0.0;

		R_Matrix[3][0] = -val_velocity[2]/val_density;
		R_Matrix[3][1] = 0.0;
		R_Matrix[3][2] = 0.0;
		R_Matrix[3][3] = 1.0/val_density;
		R_Matrix[3][4] = 0.0;

		R_Matrix[4][0] = 0.5*gm1*sqvel/val_pressure - Gamma/val_density;
		R_Matrix[4][1] = -gm1*val_velocity[0]/val_pressure;
		R_Matrix[4][2] = -gm1*val_velocity[1]/val_pressure;
		R_Matrix[4][3] = -gm1*val_velocity[2]/val_pressure;
		R_Matrix[4][4] = gm1/val_pressure;

	}

}



void CNumerics::GetRMatrix(su2double val_soundspeed, su2double val_density, su2double* val_normal, su2double **R_Matrix) {

	if (nDim == 2) {



//		R_Matrix[0][0] = 1.0;
//		R_Matrix[0][1] = 0.0;
//		R_Matrix[0][2] = 0.5*val_density/val_soundspeed;
//		R_Matrix[0][3] = 0.5*val_density/val_soundspeed;
//
//		R_Matrix[1][0] = 0.0;
//		R_Matrix[1][1] = val_normal[1];
//		R_Matrix[1][2] = 0.5*val_normal[0];
//		R_Matrix[1][3] = -0.5*val_normal[0];
//
//		R_Matrix[2][0] = 0.0;
//		R_Matrix[2][1] = -val_normal[0];
//		R_Matrix[2][2] = 0.5*val_normal[1];
//		R_Matrix[2][3] = -0.5*val_normal[1];
//
//		R_Matrix[3][0] = 0.0;
//		R_Matrix[3][1] = 0.0;
//		R_Matrix[3][2] = 0.5*val_density*val_soundspeed;
//		R_Matrix[3][3] = 0.5*val_density*val_soundspeed;

		su2double cc, rhoc;
		cc = val_soundspeed*val_soundspeed;
		rhoc = val_density/val_soundspeed;

		R_Matrix[0][0] = -1.0/cc;
		R_Matrix[0][1] = 0.0;
		R_Matrix[0][2] = 0.5/cc;
		R_Matrix[0][3] = 0.5/cc;

		R_Matrix[1][0] = 0.0;
		R_Matrix[1][1] = 0.0;
		R_Matrix[1][2] = 0.5/rhoc;
		R_Matrix[1][3] = -0.5/rhoc;

		R_Matrix[2][0] = 0.0;
		R_Matrix[2][1] = 1.0/rhoc;
		R_Matrix[2][2] = 0.0;
		R_Matrix[2][3] = 0.0;

		R_Matrix[3][0] = 0.0;
		R_Matrix[3][1] = 0.0;
		R_Matrix[3][2] = 0.5;
		R_Matrix[3][3] = 0.5;

	}
	else {

//		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
//
//		R_Matrix[0][0] =  0.5*gm1*sqvel;
//		R_Matrix[0][1] = -val_velocity[0]*gm1;
//		R_Matrix[0][2] = -val_velocity[1]*gm1;
//		R_Matrix[0][3] = -val_velocity[2]*gm1;
//		R_Matrix[0][4] = gm1;
//
//		R_Matrix[1][0] = -val_velocity[0]/val_density;
//		R_Matrix[1][1] = 1.0/val_density;
//		R_Matrix[1][2] = 0.0;
//		R_Matrix[1][3] = 0.0;
//		R_Matrix[1][4] = 0.0;
//
//		R_Matrix[2][0] = -val_velocity[1]/val_density;
//		R_Matrix[2][1] = 0.0;
//		R_Matrix[2][2] = 1.0/val_density;
//		R_Matrix[2][3] = 0.0;
//		R_Matrix[2][4] = 0.0;
//
//		R_Matrix[3][0] = -val_velocity[2]/val_density;
//		R_Matrix[3][1] = 0.0;
//		R_Matrix[3][2] = 0.0;
//		R_Matrix[3][3] = 1.0/val_density;
//		R_Matrix[3][4] = 0.0;
//
//		R_Matrix[4][0] = 0.5*gm1*sqvel/val_pressure - Gamma/val_density;
//		R_Matrix[4][1] = -gm1*val_velocity[0]/val_pressure;
//		R_Matrix[4][2] = -gm1*val_velocity[1]/val_pressure;
//		R_Matrix[4][3] = -gm1*val_velocity[2]/val_pressure;
//		R_Matrix[4][4] = gm1/val_pressure;

	}

}

void CNumerics::GetLMatrix(su2double val_soundspeed, su2double val_density, su2double* val_normal, su2double **L_Matrix) {

	if (nDim == 2) {



		L_Matrix[0][0] = 1.0;
		L_Matrix[0][1] = 0.0;
		L_Matrix[0][2] = 0.0;
		L_Matrix[0][3] = -0.5/val_soundspeed/val_soundspeed;

		L_Matrix[1][0] = 0.0;
		L_Matrix[1][1] = val_normal[1];
		L_Matrix[1][2] = -val_normal[0];
		L_Matrix[1][3] = 0.0;

		L_Matrix[2][0] = 0.0;
		L_Matrix[2][1] = val_normal[0];
		L_Matrix[2][2] = val_normal[1];
		L_Matrix[2][3] = 1.0/val_density/val_soundspeed;

		L_Matrix[3][0] = 0.0;
		L_Matrix[3][1] = -val_normal[0];
		L_Matrix[3][2] = -val_normal[1];
		L_Matrix[3][3] = 1.0/val_density/val_soundspeed;
	}
	else {

//		sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
//
//		R_Matrix[0][0] =  0.5*gm1*sqvel;
//		R_Matrix[0][1] = -val_velocity[0]*gm1;
//		R_Matrix[0][2] = -val_velocity[1]*gm1;
//		R_Matrix[0][3] = -val_velocity[2]*gm1;
//		R_Matrix[0][4] = gm1;
//
//		R_Matrix[1][0] = -val_velocity[0]/val_density;
//		R_Matrix[1][1] = 1.0/val_density;
//		R_Matrix[1][2] = 0.0;
//		R_Matrix[1][3] = 0.0;
//		R_Matrix[1][4] = 0.0;
//
//		R_Matrix[2][0] = -val_velocity[1]/val_density;
//		R_Matrix[2][1] = 0.0;
//		R_Matrix[2][2] = 1.0/val_density;
//		R_Matrix[2][3] = 0.0;
//		R_Matrix[2][4] = 0.0;
//
//		R_Matrix[3][0] = -val_velocity[2]/val_density;
//		R_Matrix[3][1] = 0.0;
//		R_Matrix[3][2] = 0.0;
//		R_Matrix[3][3] = 1.0/val_density;
//		R_Matrix[3][4] = 0.0;
//
//		R_Matrix[4][0] = 0.5*gm1*sqvel/val_pressure - Gamma/val_density;
//		R_Matrix[4][1] = -gm1*val_velocity[0]/val_pressure;
//		R_Matrix[4][2] = -gm1*val_velocity[1]/val_pressure;
//		R_Matrix[4][3] = -gm1*val_velocity[2]/val_pressure;
//		R_Matrix[4][4] = gm1/val_pressure;

	}

}

void CNumerics::GetPrecondJacobian(su2double Beta2, su2double r_hat, su2double s_hat, su2double t_hat, su2double rB2a2, su2double* Lambda, su2double *val_normal,
		su2double **val_absPeJac) {

	su2double lam1, lam2, lam3, lam4;
	lam1 = Lambda[0]; lam2 = Lambda[1]; lam3 = Lambda[2]; lam4 = Lambda[3];

	if (nDim == 2) {

		val_absPeJac[0][0] =  lam3*s_hat/(2.0*t_hat) - lam4*r_hat/(2.0*t_hat);
		val_absPeJac[0][1] = -lam3*rB2a2*val_normal[0]/(2.0*t_hat) + lam4*rB2a2*val_normal[0]/(2.0*t_hat);
		val_absPeJac[0][2] = -lam3*rB2a2*val_normal[1]/(2.0*t_hat) + lam4*rB2a2*val_normal[1]/(2.0*t_hat);
		val_absPeJac[0][3] =  0.0;

		val_absPeJac[1][0] = r_hat*val_normal[0]*lam3*s_hat/(2.0*t_hat*rB2a2) + s_hat*val_normal[0]*lam4*(-r_hat)/(2.0*t_hat*rB2a2);
		val_absPeJac[1][1] = lam2*(val_normal[1]*val_normal[1]) - lam3*r_hat*val_normal[0]*val_normal[0]/(2.0*t_hat) + lam4*s_hat*val_normal[0]*val_normal[0]/(2.0*t_hat);
		val_absPeJac[1][2] = -lam2*val_normal[0]*val_normal[1] - lam3*r_hat*val_normal[0]*val_normal[1]/(2.0*t_hat) + lam4*s_hat*val_normal[0]*val_normal[1]/(2.0*t_hat);
		val_absPeJac[1][3] = 0.0;

		val_absPeJac[2][0] = lam3*r_hat*val_normal[1]*s_hat/(2.0*t_hat*rB2a2) - s_hat*val_normal[1]*lam4*r_hat/(2.0*t_hat*rB2a2);
		val_absPeJac[2][1] = -val_normal[0]*val_normal[1]*lam2 - r_hat*val_normal[1]*val_normal[0]*lam3/(2.0*t_hat) + s_hat*val_normal[0]*val_normal[1]*lam4/(2.0*t_hat);
		val_absPeJac[2][2] = val_normal[0]*val_normal[0]*lam2 -r_hat*val_normal[1]*val_normal[1]*lam3/(2.0*t_hat) + s_hat*val_normal[1]*val_normal[1]*lam4/(2.0*t_hat);
		val_absPeJac[2][3] = 0.0;

		val_absPeJac[3][0] = 0.0;
		val_absPeJac[3][1] = 0.0;
		val_absPeJac[3][2] = 0.0;
		val_absPeJac[3][3] = lam1;

	}
	else {

		su2double lam5 = Lambda[4];

		val_absPeJac[0][0] =  lam4*s_hat/(2.0*t_hat) - lam5*r_hat/(2.0*t_hat);
		val_absPeJac[0][1] = -lam4*rB2a2*val_normal[0]/(2.0*t_hat) + lam5*rB2a2*val_normal[0]/(2.0*t_hat);
		val_absPeJac[0][2] = -lam4*rB2a2*val_normal[1]/(2.0*t_hat) + lam5*rB2a2*val_normal[1]/(2.0*t_hat);
		val_absPeJac[0][3] = -lam4*rB2a2*val_normal[2]/(2.0*t_hat) + lam5*rB2a2*val_normal[2]/(2.0*t_hat);
		val_absPeJac[0][4] =  0.0;

		val_absPeJac[1][0] = r_hat*val_normal[0]*lam4*s_hat/(2.0*t_hat*rB2a2) + s_hat*val_normal[0]*lam5*(-r_hat)/(2.0*t_hat*rB2a2);
		val_absPeJac[1][1] = lam2*(val_normal[2]*val_normal[2] + val_normal[1]*val_normal[1]) - lam4*r_hat*val_normal[0]*val_normal[0]/(2.0*t_hat) + lam5*s_hat*val_normal[0]*val_normal[0]/(2.0*t_hat);
		val_absPeJac[1][2] = -lam2*val_normal[0]*val_normal[1] - lam4*r_hat*val_normal[0]*val_normal[1]/(2.0*t_hat) + lam5*s_hat*val_normal[0]*val_normal[1]/(2.0*t_hat);
		val_absPeJac[1][3] = -lam2*val_normal[0]*val_normal[2] - lam4*r_hat*val_normal[0]*val_normal[2]/(2.0*t_hat) + lam5*s_hat*val_normal[0]*val_normal[2]/(2.0*t_hat);
		val_absPeJac[1][4] = 0.0;

		val_absPeJac[2][0] = lam4*r_hat*val_normal[1]*s_hat/(2.0*t_hat*rB2a2) - s_hat*val_normal[1]*lam5*r_hat/(2.0*t_hat*rB2a2);
		val_absPeJac[2][1] = -val_normal[0]*val_normal[1]*lam2 - r_hat*val_normal[1]*val_normal[0]*lam4/(2.0*t_hat) + s_hat*val_normal[0]*val_normal[1]*lam5/(2.0*t_hat);
		val_absPeJac[2][2] = val_normal[0]*val_normal[0]*lam2 + val_normal[2]*val_normal[2]*lam3 -r_hat*val_normal[1]*val_normal[1]*lam4/(2.0*t_hat) + s_hat*val_normal[1]*val_normal[1]*lam5/(2.0*t_hat);
		val_absPeJac[2][3] = -val_normal[2]*val_normal[1]*lam2 - r_hat*val_normal[2]*val_normal[1]*lam4/(2.0*t_hat) + s_hat*lam5*val_normal[1]*val_normal[2]/(2.0*t_hat);
		val_absPeJac[2][4] = 0.0;

		val_absPeJac[3][0] = r_hat*s_hat*val_normal[2]*lam4/(2.0*t_hat*rB2a2) - r_hat*s_hat*val_normal[2]*lam5/(2.0*t_hat*rB2a2);
		val_absPeJac[3][1] = -val_normal[0]*val_normal[2]*lam3 - lam4*val_normal[0]*val_normal[2]*r_hat/(2.0*t_hat) + lam5*val_normal[0]*val_normal[2]*s_hat/(2.0*t_hat);
		val_absPeJac[3][2] = -val_normal[1]*val_normal[2]*lam3 - lam4*val_normal[1]*val_normal[2]*r_hat/(2.0*t_hat) + lam5*val_normal[1]*val_normal[2]*s_hat/(2.0*t_hat);
		val_absPeJac[3][3] = (val_normal[1]*val_normal[1] + val_normal[0]*val_normal[0])*lam3 - lam4*val_normal[2]*val_normal[2]*r_hat/(2.0*t_hat) + lam5*val_normal[2]*val_normal[2]*s_hat/(2.0*t_hat);
		val_absPeJac[3][4] = 0.0;

		val_absPeJac[4][0] = 0.0;
		val_absPeJac[4][1] = 0.0;
		val_absPeJac[4][2] = 0.0;
		val_absPeJac[4][3] = 0.0;
		val_absPeJac[4][4] = lam1;

	}

}

void CNumerics::GetPArtCompMatrix(su2double *val_density, su2double *val_velocity, su2double *val_betainc2,
		su2double *val_normal, su2double **val_p_tensor) {
	su2double a, a2, Projvel, area2, sx, sy, sz = 0.0, u, v, w = 0.0, factor = 0.0;

	sx = val_normal[0]; sy = val_normal[1]; u = val_velocity[0]; v = val_velocity[1];
    if (nDim == 3) { sz = val_normal[2]; w = val_velocity[2]; }
	Projvel = u*sx + v*sy; area2 = sx*sx + sy*sy;
    if (nDim == 3) { Projvel += w*sz; area2 += sz*sz; }
	a2 = Projvel*Projvel + ((*val_betainc2)/(*val_density))*area2; a = sqrt(a2);
	factor = 1/(2.0*((*val_betainc2)/(*val_density))*a2);

    if (nDim == 2) {
		val_p_tensor[0][0] = 0.0;
		val_p_tensor[0][1] = factor*((*val_betainc2)/(*val_density))*a;
		val_p_tensor[0][2] = -factor*((*val_betainc2)/(*val_density))*a;
        
		val_p_tensor[1][0] = -factor*2.0*sy*((*val_betainc2)/(*val_density));
		val_p_tensor[1][1] = factor*(u*(a+Projvel) + sx*((*val_betainc2)/(*val_density)));
		val_p_tensor[1][2] = factor*(u*(Projvel-a) + sx*((*val_betainc2)/(*val_density)));
        
		val_p_tensor[2][0] = factor*2.0*sx*((*val_betainc2)/(*val_density));
		val_p_tensor[2][1] = factor*(v*(a+Projvel) + sy*((*val_betainc2)/(*val_density)));
		val_p_tensor[2][2] = factor*(v*(Projvel-a) + sy*((*val_betainc2)/(*val_density)));
	}
	else {
		val_p_tensor[0][0] = 0.0;
		val_p_tensor[0][1] = 0.0;
		val_p_tensor[0][2] = ((*val_betainc2)/(*val_density))*a;
		val_p_tensor[0][3] = -((*val_betainc2)/(*val_density))*a;

		val_p_tensor[1][0] = -sz;
		val_p_tensor[1][1] = -sy;
		val_p_tensor[1][2] = u*(Projvel+a) + sx*((*val_betainc2)/(*val_density));
		val_p_tensor[1][3] = u*(Projvel-a) + sx*((*val_betainc2)/(*val_density));

		val_p_tensor[2][0] = 0.0;
		val_p_tensor[2][1] = sx;
		val_p_tensor[2][2] = v*(Projvel+a) + sy*((*val_betainc2)/(*val_density));
		val_p_tensor[2][3] = v*(Projvel-a) + sy*((*val_betainc2)/(*val_density));

		val_p_tensor[3][0] = sx;
		val_p_tensor[3][1] = 0.0;
		val_p_tensor[3][2] = w*(Projvel+a) + sz*((*val_betainc2)/(*val_density));
		val_p_tensor[3][3] = w*(Projvel-a) + sz*((*val_betainc2)/(*val_density));
	}

}

void CNumerics::GetPArtCompMatrix_inv(su2double *val_density, su2double *val_velocity, su2double *val_betainc2,
		su2double *val_normal, su2double **val_invp_tensor) {
	su2double a, a2, Projvel, area2, sx, sy, sz = 0.0, u, v, w = 0.0;

	sx = val_normal[0]; sy = val_normal[1]; u = val_velocity[0]; v = val_velocity[1];
    if (nDim == 3) { sz = val_normal[2]; w = val_velocity[2];}
	Projvel = u*sx + v*sy; area2 = sx*sx + sy*sy;
    if (nDim == 3) { Projvel += w*sz; area2 += sz*sz; }
	a2 = Projvel*Projvel + ((*val_betainc2)/(*val_density))*area2; a = sqrt(a2);

    if (nDim == 2) {
		val_invp_tensor[0][0] = (sy*u-sx*v);
		val_invp_tensor[0][1] = -v*Projvel-sy*((*val_betainc2)/(*val_density));
		val_invp_tensor[0][2] = u*Projvel+sx*((*val_betainc2)/(*val_density));
        
		val_invp_tensor[1][0] = (a-Projvel);
		val_invp_tensor[1][1] = ((*val_betainc2)/(*val_density))*sx;
		val_invp_tensor[1][2] = ((*val_betainc2)/(*val_density))*sy;
        
		val_invp_tensor[2][0] = (-a-Projvel);
		val_invp_tensor[2][1] = ((*val_betainc2)/(*val_density))*sx;
		val_invp_tensor[2][2] = ((*val_betainc2)/(*val_density))*sy;
	}
	else {
		val_invp_tensor[0][0] = (sz*Projvel-area2*w)/(sx*a2);
		val_invp_tensor[0][1] = -(w*Projvel+sz*((*val_betainc2)/(*val_density)))/a2;
		val_invp_tensor[0][2] = -sy*(w*Projvel+sz*((*val_betainc2)/(*val_density)))/(sx*a2);
		val_invp_tensor[0][3] = ((sx*u+sy*v)*Projvel+(sx*sx+sy*sy)*((*val_betainc2)/(*val_density)))/(sx*a2);

		val_invp_tensor[1][0] = (sy*Projvel-area2*v)/(sx*a2);
		val_invp_tensor[1][1] = -(v*Projvel+sy*((*val_betainc2)/(*val_density)))/a2;
		val_invp_tensor[1][2] = ((sx*u+sz*w)*Projvel+(sx*sx+sz*sz)*((*val_betainc2)/(*val_density)))/(sx*a2);
		val_invp_tensor[1][3] = -sz*(v*Projvel+sy*((*val_betainc2)/(*val_density)))/(sx*a2);

		val_invp_tensor[2][0] = -(Projvel-a)/(2.0*a2*((*val_betainc2)/(*val_density)));
		val_invp_tensor[2][1] = sx/(2.0*a2);
		val_invp_tensor[2][2] = sy/(2.0*a2);
		val_invp_tensor[2][3] = sz/(2.0*a2);

		val_invp_tensor[3][0] = -(Projvel+a)/(2.0*a2*((*val_betainc2)/(*val_density)));
		val_invp_tensor[3][1] = sx/(2.0*a2);
		val_invp_tensor[3][2] = sy/(2.0*a2);
		val_invp_tensor[3][3] = sz/(2.0*a2);
	}

}

void CNumerics::GetPArtComp_FreeSurf_Matrix(su2double *val_density, su2double *val_ddensity, su2double *val_velocity, su2double *val_betainc2, su2double *val_levelset, su2double *val_normal, su2double **val_p_tensor) {
  
	su2double a = 0.0, b = 0.0, c = 0.0, d = 0.0, area2 = 0.0, e2 = 0.0, f = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, u = 0.0, v = 0.0, w = 0.0;
  
  a = (*val_betainc2)/(*val_density);
  b = (*val_levelset)/(*val_density);
  c = (*val_ddensity);
  
  if (nDim == 2) {
    
    nx = val_normal[0];   ny = val_normal[1];   area2 = nx*nx + ny*ny;
    u = val_velocity[0];  v = val_velocity[1];  d = u*nx + v*ny;
    e2 = (2.0*d - b*c*d)*(2.0*d - b*c*d);
    f = sqrt(4.0*a*area2 + e2);
    
    val_p_tensor[0][0] = 0;
    val_p_tensor[0][1] = 0;
    val_p_tensor[0][2] = (d*d*(1.0 - b*c) + 2.0*a*area2 + d*d + d*f)/(2.0*b*area2);
    val_p_tensor[0][3] = (d*d*(1.0 - b*c) + 2.0*a*area2 + d*d - d*f)/(2.0*b*area2);
    
    val_p_tensor[1][0] = (c*d)/nx;
    val_p_tensor[1][1] = -(ny/nx);
    val_p_tensor[1][2] = (d*nx*(b*c - 1.0) + nx*nx*u + 2.0*ny*ny*u - nx*ny*v - nx*f)/(2*b*area2);
    val_p_tensor[1][3] = (d*nx*(b*c - 1.0) + nx*nx*u + 2.0*ny*ny*u - nx*ny*v + nx*f)/(2*b*area2);
    
    val_p_tensor[2][0] = 0.0;
    val_p_tensor[2][1] = 1.0;
    val_p_tensor[2][2] = (d*ny*(b*c - 1.0) - nx*ny*u + 2.0*nx*nx*v + ny*ny*v - ny*f)/(2*b*area2);
    val_p_tensor[2][3] = (d*ny*(b*c - 1.0) - nx*ny*u + 2.0*nx*nx*v + ny*ny*v + ny*f)/(2*b*area2);
    
    val_p_tensor[3][0] = 1.0;
    val_p_tensor[3][1] = 0.0;
    val_p_tensor[3][2] = 1.0;
    val_p_tensor[3][3] = 1.0;
    
	}
	else {
  
    nx = val_normal[0];   ny = val_normal[1];   nz = val_normal[2];   area2 = nx*nx + ny*ny + nz*nz;
    u = val_velocity[0];  v = val_velocity[1];  w = val_velocity[2];  d = u*nx + v*ny + w*nz;
    e2 = (2.0*d - b*c*d)*(2.0*d - b*c*d);
    f = sqrt(4.0*a*area2 + e2);
  
    val_p_tensor[0][0] = 0.0;
    val_p_tensor[0][1] = 0.0;
    val_p_tensor[0][2] = 0.0;
    val_p_tensor[0][3] = -((a*(b*c*d + f))/(b*(2.0*d - b*c*d - f)));
    val_p_tensor[0][4] = (a*(- b*c*d + f))/(b*(2.0*d - b*c*d + f));
      
    val_p_tensor[1][0] = (c*d)/nx;
    val_p_tensor[1][1] = -(nz/nx);
    val_p_tensor[1][2] = -(ny/nx);
    val_p_tensor[1][3] = -((-2.0*a*nx + b*c*d*u - 2.0*u*d + u*f)/(b*(2.0*d - b*c*d - f)));
    val_p_tensor[1][4] = -((-2.0*a*nx + b*c*d*u - 2.0*u*d - u*f)/(b*(2.0*d - b*c*d + f)));
      
    val_p_tensor[2][0] = 0.0;
    val_p_tensor[2][1] = 0.0;
    val_p_tensor[2][2] = 1.0;
    val_p_tensor[2][3] = -((-2.0*a*ny + b*c*d*v - 2.0*v*d + v*f)/(b*(2.0*d - b*c*d - f)));
    val_p_tensor[2][4] = -((-2.0*a*ny + b*c*d*v - 2.0*v*d - v*f)/(b*(2.0*d - b*c*d + f)));
      
    val_p_tensor[3][0] = 0.0;
    val_p_tensor[3][1] = 1.0;
    val_p_tensor[3][2] = 0.0;
    val_p_tensor[3][3] = -((-2.0*a*nz + b*c*d*w - 2.0*w*d + w*f)/(b*(2.0*d - b*c*d - f)));
    val_p_tensor[3][4] = -((-2.0*a*nz + b*c*d*w - 2.0*w*d - w*f)/(b*(2.0*d - b*c*d + f)));

    val_p_tensor[4][0] = 1.0;
    val_p_tensor[4][1] = 0.0;
    val_p_tensor[4][2] = 0.0;
    val_p_tensor[4][3] = 1.0;
    val_p_tensor[4][4] = 1.0;
    
  }
  
}

void CNumerics::GetPArtComp_FreeSurf_Matrix_inv(su2double *val_density, su2double *val_ddensity, su2double *val_velocity, su2double *val_betainc2, su2double *val_levelset,
                                      su2double *val_normal, su2double **val_invp_tensor) {
  
	su2double a = 0.0, b = 0.0, c = 0.0, d = 0.0, area2 = 0.0, e2 = 0.0, f = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, u = 0.0, v = 0.0, w = 0.0;
  
  a = (*val_betainc2)/(*val_density);
  b = (*val_levelset)/(*val_density);
  c = (*val_ddensity);
  
  if (nDim == 2) {
    
    nx = val_normal[0];   ny = val_normal[1];   area2 = nx*nx + ny*ny;
    u = val_velocity[0];  v = val_velocity[1];  d = u*nx + v*ny;
    e2 = (2.0*d - b*c*d)*(2.0*d - b*c*d);
    f = sqrt(4.0*a*area2 + e2);

    val_invp_tensor[0][0] = -((b*area2)/(a*area2 + d*d*(1.0 - b*c)));
    val_invp_tensor[0][1] = -((b*d*nx)/(a*area2 + d*d*(1.0 - b*c)));
    val_invp_tensor[0][2] = -((b*d*ny)/(a*area2 + d*d*(1.0 - b*c)));
    val_invp_tensor[0][3] = ( a*area2 + d*d)/(a*area2 + d*d*(1.0 - b*c));
    
    val_invp_tensor[1][0] = (-b*c*d*ny + nx*(ny*u - nx*v))/(a*area2 + d*d*(1.0-b*c));
    val_invp_tensor[1][1] = -((nx*(a*ny + d*v))/(a*area2 + d*d*(1.0-b*c)));
    val_invp_tensor[1][2] = (-b*c*d*d + nx*(a*nx + d*u))/(a*area2 + d*d*(1.0-b*c));
    val_invp_tensor[1][3] = (c*d*(a*ny + d*v))/(a*area2 + d*d*(1.0-b*c));
    
    val_invp_tensor[2][0] = (b*area2*(-b*c*d + f))/(2.0*(-b*c*d*d + a*area2 + d*d)*f);
    val_invp_tensor[2][1] = -((b*nx*(-(2.0 - b*c)*d*d - 2.0*a*area2 + d*f))/(2.0*(- a*area2 - d*d*(1.0-b*c))*f));
    val_invp_tensor[2][2] = -((b*ny*(-(2.0 - b*c)*d*d - 2.0*a*area2 + d*f))/(2.0*(- a*area2 - d*d*(1.0-b*c))*f));
    val_invp_tensor[2][3] = (b*c*d*(-(2.0 - b*c)*d*d - 2.0*a*area2 + d*f))/(2.0*(- a*area2 - d*d*(1.0-b*c))*f);
    
    val_invp_tensor[3][0] = (b*area2*(b*c*d + f))/(2.0*(-b*c*d*d + a*area2 + d*d)*f);
    val_invp_tensor[3][1] = -((b*nx*((2.0 - b*c)*d*d + 2.0*a*area2 + d*f))/(2.0*(- a*area2 - d*d*(1.0-b*c))*f));
    val_invp_tensor[3][2] = -((b*ny*((2.0 - b*c)*d*d + 2.0*a*area2 + d*f))/(2.0*(- a*area2 - d*d*(1.0-b*c))*f));
    val_invp_tensor[3][3] = (b*c*d*((2.0 - b*c)*d*d + 2.0*a*area2 + d*f))/(2.0*(- a*area2 - d*d*(1.0-b*c))*f);
    
	}
	else {
  
    nx = val_normal[0];   ny = val_normal[1];   nz = val_normal[2];   area2 = nx*nx + ny*ny + nz*nz;
    u = val_velocity[0];  v = val_velocity[1];  w = val_velocity[2];  d = u*nx + v*ny + w*nz;
    e2 = (2.0*d - b*c*d)*(2.0*d - b*c*d);
    f = sqrt(4.0*a*area2 + e2);
    
    val_invp_tensor[0][0] = (b*area2)/(b*c*d*d - a*area2 - d*d);
    val_invp_tensor[0][1] = -((b*d*nx)/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[0][2] = -((b*d*ny)/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[0][3] = -((b*d*nz)/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[0][4] = (a*area2 + d*d)/(-b*c*d*d +a*area2 + d*d);
    
    val_invp_tensor[1][0] = (-b*c*d*nz +nx*nz*u - nx*nx*w + ny*(nz*v - ny*w))/(-b*c*d*d +a*area2 + d*d);
    val_invp_tensor[1][1] = -((nx*(a*nz + d*w))/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[1][2] = -((ny*(a*nz + d*w))/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[1][3] = (-b*c*d*d + a*(nx*nx + ny*ny) + d*(nx*u + ny*v))/(-b*c*d*d + a*area2 + d*d);
    val_invp_tensor[1][4] = (c*d*(a*nz + d*w))/(-b*c*d*d + a*area2 + d*d);
        
    val_invp_tensor[2][0] = (-b*c*d*ny + nx*ny*u - nx*nx*v + nz*(-nz*v + ny*w))/(-b*c*d*d + a*area2 + d*d);
    val_invp_tensor[2][1] =  -((nx*(a*ny + d*v))/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[2][2] = (-b*c*d*d + a*(nx*nx + nz*nz) + d*(nx*u + nz*w))/(-b*c*d*d + a*area2 + d*d);
    val_invp_tensor[2][3] = -((nz*(a*ny + d*v))/(-b*c*d*d + a*area2 + d*d));
    val_invp_tensor[2][4] = (c*d*(a*ny + d*v))/(-b*c*d*d + a*area2 + d*(nx*u + ny*v + nz*w));
    
    val_invp_tensor[3][0] = -(b*(-d + b*c*d - d + f)*(b*b*c*c*d*d + 2.0*a*area2 + d*d - 3.0*b*c*d*d + 2.0*nx*ny*u*v + 2.0*nx*nz*u*w + 2.0*ny*nz*v*w + nx*nx*u*u + ny*ny*v*v + nz*nz*w*w - b*c*d*f + d*f))/(4*a*(b*c*d*d - a*area2 - d*d)*f);
    val_invp_tensor[3][1] = (b*nx*(-d + b*c*d - d + f)*(- b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f);
    val_invp_tensor[3][2] = (b*ny*(-d + b*c*d - d + f)*(- b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f);
    val_invp_tensor[3][3] = (b*nz*(-d + b*c*d - d + f)*(- b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f);
    val_invp_tensor[3][4] = -((b*c*d*(-d + b*c*d - d + f)*(- b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f));
        
    val_invp_tensor[4][0] = -(b*(2.0*d - b*c*d + f)*(b*b*c*c*d*d + 2.0*a*area2 + d*d - 3.0*b*c*d*d + 2.0*nx*ny*u*v + 2.0*nx*nz*u*w + 2.0*ny*nz*v*w + nx*nx*u*u + ny*ny*v*v + nz*nz*w*w + b*c*d*f - d*f))/(4*a*(b*c*d*d - a*area2 - d*d)*f);
    val_invp_tensor[4][1] = -((b*nx*(b*c*d + f)*(2.0*d - b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f));
    val_invp_tensor[4][2] = -((b*ny*(b*c*d + f)*(2.0*d - b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f));
    val_invp_tensor[4][3] = -((b*nz*(b*c*d + f)*(2.0*d - b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f));
    val_invp_tensor[4][4] = (b*c*d*(b*c*d + f)*(2.0*d - b*c*d + f))/(4.0*(b*c*d*d - a*area2 - d*d)*f);

  }
  
}

void CNumerics::GetJacInviscidLambda_fabs(su2double *val_velocity, su2double val_soundspeed,
		su2double *val_normal, su2double *val_Lambda_Vector) {
	su2double ProjVelocity = 0;

	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		ProjVelocity += val_velocity[iDim]*val_normal[iDim];

	if (nDim == 3) {
		val_Lambda_Vector[0] = fabs(ProjVelocity);
		val_Lambda_Vector[1] = fabs(ProjVelocity);
		val_Lambda_Vector[2] = fabs(ProjVelocity);
		val_Lambda_Vector[3] = fabs(ProjVelocity + val_soundspeed);
		val_Lambda_Vector[4] = fabs(ProjVelocity - val_soundspeed);
	}

	if (nDim == 2) {
		val_Lambda_Vector[0] = fabs(ProjVelocity);
		val_Lambda_Vector[1] = fabs(ProjVelocity);
		val_Lambda_Vector[2] = fabs(ProjVelocity + val_soundspeed);
		val_Lambda_Vector[3] = fabs(ProjVelocity - val_soundspeed);
	}
}

void CNumerics::GetAdjViscousFlux_Jac(su2double Pressure_i, su2double Pressure_j, su2double Density_i, su2double Density_j,
                                      su2double ViscDens_i, su2double ViscDens_j, su2double *Velocity_i, su2double *Velocity_j,
                                      su2double sq_vel_i, su2double sq_vel_j,
                                      su2double XiDens_i, su2double XiDens_j, su2double **Mean_GradPhi, su2double *Mean_GradPsiE,
                                      su2double dPhiE_dn, su2double *Normal, su2double *Edge_Vector, su2double dist_ij_2, su2double *val_residual_i, su2double *val_residual_j,
                                      su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji,
                                      su2double **val_Jacobian_jj, bool implicit) {
  
  su2double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  Sigma_yz5, Sigma_5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  su2double dSigmaxx_phi1, dSigmayy_phi1, dSigmazz_phi1, dSigmaxy_phi1, dSigmaxz_phi1, dSigmayz_phi1;
  su2double dSigmaxx_phi2, dSigmayy_phi2, dSigmazz_phi2, dSigmaxy_phi2, dSigmaxz_phi2, dSigmayz_phi2;
  su2double dSigmaxx_phi3, dSigmayy_phi3, dSigmazz_phi3, dSigmaxy_phi3, dSigmaxz_phi3, dSigmayz_phi3;
  su2double dSigma5_psi5;
  unsigned short iVar, jVar;
  
  if (nDim == 3) {
    
    /*--- Residual at iPoint ---*/
    
    Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_zz = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
    Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xz = ViscDens_i * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
    Sigma_yz = ViscDens_i * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
    Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
    Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
    Sigma_zz5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_i[2] * Mean_GradPsiE[2]);
    Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
    Sigma_xz5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[0]);
    Sigma_yz5 = ViscDens_i * (Velocity_i[1] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[1]);
    Sigma_5   = XiDens_i * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
    eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
    
    val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy + Velocity_i[2] * Normal[2] * eta_zz
                           + (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
                           + (Velocity_i[0] * Normal[2] + Velocity_i[2] * Normal[0]) * eta_xz
                           + (Velocity_i[2] * Normal[1] + Velocity_i[1] * Normal[2]) * eta_yz
                           - (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
    
    val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_i[0] * Sigma_5);
    val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_i[1] * Sigma_5);
    val_residual_i[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_i[2] * Sigma_5);
    val_residual_i[4] = (Sigma_5);
    
    /*--- Computation of the Jacobians at Point i---*/
    
    if (implicit) {
      
      dSigmaxx_phi1 = -FOUR3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxx_phi3 =   TWO3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmayy_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = -FOUR3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi3 =   TWO3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmazz_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmazz_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmazz_phi3 = -FOUR3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmaxy_phi1 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmaxy_phi3 = 0;
      dSigmaxz_phi1 = -ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmaxz_phi2 = 0;
      dSigmaxz_phi3 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmayz_phi1 = 0;
      dSigmayz_phi2 = -ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmayz_phi3 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
      
      //      dSigmaxx5_psi5 = -ViscDens_i * ( FOUR3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] -  TWO3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmayy5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] + FOUR3*Velocity_i[1]*Edge_Vector[1] -  TWO3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmazz5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] + FOUR3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmaxy5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[1] + Velocity_i[1]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmaxz5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[2] + Velocity_i[2]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmayz5_psi5 = -ViscDens_i * ( Velocity_i[1]*Edge_Vector[2] + Velocity_i[2]*Edge_Vector[1] )/dist_ij_2;
      dSigma5_psi5   = -XiDens_i * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;
      
      val_Jacobian_ii[0][0] = 0;
      val_Jacobian_ii[0][1] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi1 + Velocity_i[1]*Normal[1]*dSigmayy_phi1 + Velocity_i[2]*Normal[2]*dSigmazz_phi1
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi1
                                + (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi1
                                + (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi1 );
      val_Jacobian_ii[0][2] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi2 + Velocity_i[1]*Normal[1]*dSigmayy_phi2 + Velocity_i[2]*Normal[2]*dSigmazz_phi2
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi2
                                + (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi2
                                + (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi2 );
      val_Jacobian_ii[0][3] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi3 + Velocity_i[1]*Normal[1]*dSigmayy_phi3 + Velocity_i[2]*Normal[2]*dSigmazz_phi3
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi3
                                + (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi3
                                + (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi3 );
      val_Jacobian_ii[0][4] = (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * dSigma5_psi5;
      
      val_Jacobian_ii[1][0] = 0;
      val_Jacobian_ii[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1 + Normal[2]*dSigmaxz_phi1;
      val_Jacobian_ii[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2 + Normal[2]*dSigmaxz_phi2;
      val_Jacobian_ii[1][3] = Normal[0]*dSigmaxx_phi3 + Normal[1]*dSigmaxy_phi3 + Normal[2]*dSigmaxz_phi3;
      val_Jacobian_ii[1][4] = -Velocity_i[0]*dSigma5_psi5;
      
      val_Jacobian_ii[2][0] = 0;
      val_Jacobian_ii[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1 + Normal[2]*dSigmayz_phi1;
      val_Jacobian_ii[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2 + Normal[2]*dSigmayz_phi2;
      val_Jacobian_ii[2][3] = Normal[0]*dSigmaxy_phi3 + Normal[1]*dSigmayy_phi3 + Normal[2]*dSigmayz_phi3;
      val_Jacobian_ii[2][4] = -Velocity_i[1]*dSigma5_psi5;
      
      val_Jacobian_ii[3][0] = 0;
      val_Jacobian_ii[3][1] = Normal[0]*dSigmaxz_phi1 + Normal[1]*dSigmayz_phi1 + Normal[2]*dSigmazz_phi1;
      val_Jacobian_ii[3][2] = Normal[0]*dSigmaxz_phi2 + Normal[1]*dSigmayz_phi2 + Normal[2]*dSigmazz_phi2;
      val_Jacobian_ii[3][3] = Normal[0]*dSigmaxz_phi3 + Normal[1]*dSigmayz_phi3 + Normal[2]*dSigmazz_phi3;
      val_Jacobian_ii[3][4] = -Velocity_i[2]*dSigma5_psi5;
      
      val_Jacobian_ii[4][0] = 0;
      val_Jacobian_ii[4][1] = 0;
      val_Jacobian_ii[4][2] = 0;
      val_Jacobian_ii[4][3] = 0;
      val_Jacobian_ii[4][4] = dSigma5_psi5;
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ij[iVar][jVar] = -val_Jacobian_ii[iVar][jVar];
    }
    
    /*--- Residual at jPoint ---*/
    
    Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_zz = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
    Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xz = ViscDens_j * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
    Sigma_yz = ViscDens_j * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
    Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
    Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
    Sigma_zz5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_j[2] * Mean_GradPsiE[2]);
    Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
    Sigma_xz5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[0]);
    Sigma_yz5 = ViscDens_j * (Velocity_j[1] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[1]);
    Sigma_5   = XiDens_j * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
    eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
    
    val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy + Velocity_j[2] * Normal[2] * eta_zz
                           + (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
                           + (Velocity_j[0] * Normal[2] + Velocity_j[2] * Normal[0]) * eta_xz
                           + (Velocity_j[2] * Normal[1] + Velocity_j[1] * Normal[2]) * eta_yz
                           - (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
    val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_j[0] * Sigma_5);
    val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_j[1] * Sigma_5);
    val_residual_j[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_j[2] * Sigma_5);
    val_residual_j[4] = (Sigma_5);
    
    /*--- Computation of the Jacobians at Point j---*/
    
    if (implicit) {
      
      dSigmaxx_phi1 = FOUR3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxx_phi3 = -TWO3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmayy_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = FOUR3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi3 = -TWO3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmazz_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmazz_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmazz_phi3 = FOUR3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmaxy_phi1 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmaxy_phi3 = 0;
      dSigmaxz_phi1 = ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmaxz_phi2 = 0;
      dSigmaxz_phi3 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmayz_phi1 = 0;
      dSigmayz_phi2 = ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmayz_phi3 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
      
      //      dSigmaxx5_psi5 = ViscDens_j * ( FOUR3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] -  TWO3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmayy5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] + FOUR3*Velocity_j[1]*Edge_Vector[1] -  TWO3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmazz5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] + FOUR3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmaxy5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[1] + Velocity_j[1]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmaxz5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[2] + Velocity_j[2]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmayz5_psi5 = ViscDens_j * ( Velocity_j[1]*Edge_Vector[2] + Velocity_j[2]*Edge_Vector[1] )/dist_ij_2;
      dSigma5_psi5   = XiDens_j * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;
      
      val_Jacobian_jj[0][0] = 0;
      val_Jacobian_jj[0][1] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi1 + Velocity_j[1]*Normal[1]*dSigmayy_phi1 + Velocity_j[2]*Normal[2]*dSigmazz_phi1
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi1
                                + (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi1
                                + (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi1 );
      val_Jacobian_jj[0][2] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi2 + Velocity_j[1]*Normal[1]*dSigmayy_phi2 + Velocity_j[2]*Normal[2]*dSigmazz_phi2
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi2
                                + (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi2
                                + (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi2 );
      val_Jacobian_jj[0][3] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi3 + Velocity_j[1]*Normal[1]*dSigmayy_phi3 + Velocity_j[2]*Normal[2]*dSigmazz_phi3
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi3
                                + (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi3
                                + (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi3 );
      val_Jacobian_jj[0][4] = (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * dSigma5_psi5;
      
      val_Jacobian_jj[1][0] = 0;
      val_Jacobian_jj[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1 + Normal[2]*dSigmaxz_phi1;
      val_Jacobian_jj[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2 + Normal[2]*dSigmaxz_phi2;
      val_Jacobian_jj[1][3] = Normal[0]*dSigmaxx_phi3 + Normal[1]*dSigmaxy_phi3 + Normal[2]*dSigmaxz_phi3;
      val_Jacobian_jj[1][4] = -Velocity_j[0]*dSigma5_psi5;
      
      val_Jacobian_jj[2][0] = 0;
      val_Jacobian_jj[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1 + Normal[2]*dSigmayz_phi1;
      val_Jacobian_jj[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2 + Normal[2]*dSigmayz_phi2;
      val_Jacobian_jj[2][3] = Normal[0]*dSigmaxy_phi3 + Normal[1]*dSigmayy_phi3 + Normal[2]*dSigmayz_phi3;
      val_Jacobian_jj[2][4] = -Velocity_j[1]*dSigma5_psi5;
      
      val_Jacobian_jj[3][0] = 0;
      val_Jacobian_jj[3][1] = Normal[0]*dSigmaxz_phi1 + Normal[1]*dSigmayz_phi1 + Normal[2]*dSigmazz_phi1;
      val_Jacobian_jj[3][2] = Normal[0]*dSigmaxz_phi2 + Normal[1]*dSigmayz_phi2 + Normal[2]*dSigmazz_phi2;
      val_Jacobian_jj[3][3] = Normal[0]*dSigmaxz_phi3 + Normal[1]*dSigmayz_phi3 + Normal[2]*dSigmazz_phi3;
      val_Jacobian_jj[3][4] = -Velocity_j[2]*dSigma5_psi5;
      
      val_Jacobian_jj[4][0] = 0;
      val_Jacobian_jj[4][1] = 0;
      val_Jacobian_jj[4][2] = 0;
      val_Jacobian_jj[4][3] = 0;
      val_Jacobian_jj[4][4] = dSigma5_psi5;
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
    }
    
  }
  else if (nDim == 2) {
    
    /*--- Residual at iPoint ---*/
    
    Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
    Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
    Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1]);
    Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1]);
    Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
    Sigma_5   = XiDens_i * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;
    
    val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy
                           + (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
                           - (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
    val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] - Velocity_i[0] * Sigma_5);
    val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] - Velocity_i[1] * Sigma_5);
    val_residual_i[3] = (Sigma_5);
    
    /*--- Computation of the Jacobians at Point i---*/
    
    if (implicit) {
      
      dSigmaxx_phi1 = -FOUR3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = -FOUR3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi1 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
      
      //      dSigmaxx5_psi5 = -ViscDens_i * ( FOUR3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmayy5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] + FOUR3*Velocity_i[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmaxy5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[1] + Velocity_i[1]*Edge_Vector[0] )/dist_ij_2;
      dSigma5_psi5   = -XiDens_i * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;
      
      val_Jacobian_ii[0][0] = 0;
      
      val_Jacobian_ii[0][1] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi1 + Velocity_i[1]*Normal[1]*dSigmayy_phi1
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi1 );
      val_Jacobian_ii[0][2] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi2 + Velocity_i[1]*Normal[1]*dSigmayy_phi2
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi2 );
      val_Jacobian_ii[0][3] = (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * dSigma5_psi5;
      
      val_Jacobian_ii[1][0] = 0;
      val_Jacobian_ii[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1;
      val_Jacobian_ii[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2;
      val_Jacobian_ii[1][3] = -Velocity_i[0]*dSigma5_psi5;
      
      val_Jacobian_ii[2][0] = 0;
      val_Jacobian_ii[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1;
      val_Jacobian_ii[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2;
      val_Jacobian_ii[2][3] = -Velocity_i[1]*dSigma5_psi5;
      
      val_Jacobian_ii[3][0] = 0;
      val_Jacobian_ii[3][1] = 0;
      val_Jacobian_ii[3][2] = 0;
      val_Jacobian_ii[3][3] = dSigma5_psi5;
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ij[iVar][jVar] = -val_Jacobian_ii[iVar][jVar];
    }
    
    /*--- Residual at jPoint ---*/
    Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
    Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
    Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1]);
    Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1]);
    Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
    Sigma_5   = XiDens_j * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;
    
    val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy
                           + (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
                           - (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
    val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1]  - Velocity_j[0] * Sigma_5);
    val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1]  - Velocity_j[1] * Sigma_5);
    val_residual_j[3] = (Sigma_5);
    
    /*--- Computation of the Jacobians at Point j---*/
    if (implicit) {
      dSigmaxx_phi1 = FOUR3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = FOUR3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi1 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
      
      //      dSigmaxx5_psi5 = ViscDens_j * ( FOUR3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmayy5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] + FOUR3*Velocity_j[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmaxy5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[1] + Velocity_j[1]*Edge_Vector[0] )/dist_ij_2;
      dSigma5_psi5   = XiDens_j * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;
      
      val_Jacobian_jj[0][0] = 0;
      val_Jacobian_jj[0][1] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi1 + Velocity_j[1]*Normal[1]*dSigmayy_phi1
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi1 );
      val_Jacobian_jj[0][2] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi2 + Velocity_j[1]*Normal[1]*dSigmayy_phi2
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi2 );
      val_Jacobian_jj[0][3] = (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * dSigma5_psi5;
      
      val_Jacobian_jj[1][0] = 0;
      val_Jacobian_jj[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1;
      val_Jacobian_jj[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2;
      val_Jacobian_jj[1][3] = -Velocity_j[0]*dSigma5_psi5;
      
      val_Jacobian_jj[2][0] = 0;
      val_Jacobian_jj[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1;
      val_Jacobian_jj[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2;
      val_Jacobian_jj[2][3] = -Velocity_j[1]*dSigma5_psi5;
      
      val_Jacobian_jj[3][0] = 0;
      val_Jacobian_jj[3][1] = 0;
      val_Jacobian_jj[3][2] = 0;
      val_Jacobian_jj[3][3] = dSigma5_psi5;
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
    }
  }
 
}

void CNumerics::GetViscousFlux(su2double *val_primvar, su2double **val_gradprimvar,
		su2double val_laminar_viscosity, su2double val_eddy_viscosity, su2double val_mach_inf) {

	su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
	su2double heat_flux_factor = Cp * (val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb);

	su2double div_vel = 0.0;
	for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
		div_vel += val_gradprimvar[iDim+1][iDim];

	for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
		for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
			tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] +
                                          val_gradprimvar[iDim+1][jDim] )
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



void CNumerics::GetViscousProjFlux(su2double *val_primvar,
									su2double **val_gradprimvar, su2double val_turb_ke,
									su2double *val_normal,
									su2double val_laminar_viscosity,
									su2double val_eddy_viscosity) {


	unsigned short iVar, iDim, jDim;
	su2double total_viscosity, heat_flux_factor, div_vel, Cp, Density;

	Density = val_primvar[nDim+2];
	total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
	heat_flux_factor = Cp * (val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb);

	div_vel = 0.0;
	for (iDim = 0 ; iDim < nDim; iDim++)
		div_vel += val_gradprimvar[iDim+1][iDim];
	for (iDim = 0 ; iDim < nDim; iDim++)
		for (jDim = 0 ; jDim < nDim; jDim++)
			tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
			- TWO3*total_viscosity*div_vel*delta[iDim][jDim]
			                                           - TWO3*Density*val_turb_ke*delta[iDim][jDim];
	/*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
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
	} else {
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
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Flux_Tensor[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
	}
}

void CNumerics::GetViscousProjFlux(su2double *val_primvar,
                                   su2double **val_gradprimvar, su2double val_turb_ke,
                                   su2double *val_normal,
                                   su2double val_laminar_viscosity,
                                   su2double val_eddy_viscosity,
                                   su2double val_thermal_conductivity,
                                   su2double val_heat_capacity_cp) {

	unsigned short iVar, iDim, jDim;
	su2double total_viscosity, heat_flux_factor, div_vel, Density;
	Density = val_primvar[nDim+2];

	total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	heat_flux_factor = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;

	div_vel = 0.0;
	for (iDim = 0 ; iDim < nDim; iDim++)
		div_vel += val_gradprimvar[iDim+1][iDim];

	for (iDim = 0 ; iDim < nDim; iDim++)
		for (jDim = 0 ; jDim < nDim; jDim++)
			tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
			- TWO3*total_viscosity*div_vel*delta[iDim][jDim]
			                                           - TWO3*Density*val_turb_ke*delta[iDim][jDim];


	/*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
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
	} else {
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

	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Flux_Tensor[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
	}

}

void CNumerics::GetViscousArtCompProjFlux(su2double **val_gradprimvar, su2double *val_normal, su2double val_laminar_viscosity,
		su2double val_eddy_viscosity) {
	unsigned short iVar, iDim;
	su2double total_viscosity;
	
	total_viscosity = (val_laminar_viscosity + val_eddy_viscosity);

	if (nDim == 3) {
		Flux_Tensor[0][0] = 0.0;
		Flux_Tensor[1][0] = total_viscosity * val_gradprimvar[1][0];
		Flux_Tensor[2][0] = total_viscosity * val_gradprimvar[2][0];
		Flux_Tensor[3][0] = total_viscosity * val_gradprimvar[3][0];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = total_viscosity * val_gradprimvar[1][1];
		Flux_Tensor[2][1] = total_viscosity * val_gradprimvar[2][1];
		Flux_Tensor[3][1] = total_viscosity * val_gradprimvar[3][1];

		Flux_Tensor[0][2] = 0.0;
		Flux_Tensor[1][2] = total_viscosity * val_gradprimvar[1][2];
		Flux_Tensor[2][2] = total_viscosity * val_gradprimvar[2][2];
		Flux_Tensor[3][2] = total_viscosity * val_gradprimvar[3][2];
	}

	if (nDim == 2) {
		Flux_Tensor[0][0] = 0.0;
		Flux_Tensor[1][0] = total_viscosity * val_gradprimvar[1][0];
		Flux_Tensor[2][0] = total_viscosity * val_gradprimvar[2][0];

		Flux_Tensor[0][1] = 0.0;
		Flux_Tensor[1][1] = total_viscosity * val_gradprimvar[1][1];
		Flux_Tensor[2][1] = total_viscosity * val_gradprimvar[2][1];
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Flux_Tensor[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
	}
}

void CNumerics::GetViscousProjJacs(su2double *val_Mean_PrimVar, su2double val_laminar_viscosity,
		su2double val_eddy_viscosity, su2double val_dist_ij, su2double *val_normal, su2double val_dS,
		su2double *val_Proj_Visc_Flux, su2double **val_Proj_Jac_Tensor_i, su2double **val_Proj_Jac_Tensor_j) {
	unsigned short iDim, iVar, jVar;

	su2double theta = 0.0, sqvel = 0.0, proj_viscousflux_vel = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim++) {
		theta += val_normal[iDim]*val_normal[iDim];
		sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
		proj_viscousflux_vel += val_Proj_Visc_Flux[iDim+1]*val_Mean_PrimVar[iDim+1];
	}
  
	su2double phi = 0.5*(Gamma-1.0)*sqvel;
	su2double Density = val_Mean_PrimVar[nDim+2];
	su2double Pressure = val_Mean_PrimVar[nDim+1];
	su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	su2double heat_flux_factor = val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb;
	su2double cpoR = Gamma/(Gamma-1.0); // cp over R
	su2double factor = total_viscosity*val_dS/(Density*val_dist_ij);
	su2double phi_rho = -cpoR*heat_flux_factor*Pressure/(Density*Density);
	su2double phi_p = cpoR*heat_flux_factor/(Density);
	su2double rhoovisc = Density/(total_viscosity); // rho over viscosity

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      val_Proj_Jac_Tensor_i[iVar][jVar] = 0.0;
      val_Proj_Jac_Tensor_j[iVar][jVar] = 0.0;
    }
  }
  
	if (nDim == 2) {
    
		su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;

		su2double etaz = val_normal[0]*val_normal[1]/3.0;

		su2double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz;
		su2double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay;

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

		val_Proj_Jac_Tensor_i[3][0] = -factor*(rhoovisc*theta*(phi_rho+phi*phi_p)- (pix*val_Mean_PrimVar[1]+piy*val_Mean_PrimVar[2]));
		val_Proj_Jac_Tensor_i[3][1] = -factor*(pix-rhoovisc*theta*phi_p*(Gamma-1.0)*val_Mean_PrimVar[1]);
		val_Proj_Jac_Tensor_i[3][2] = -factor*(piy-rhoovisc*theta*phi_p*(Gamma-1.0)*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[3][3] = -factor*((Gamma-1.0)*rhoovisc*theta*phi_p);
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		factor = 0.5/Density;
		val_Proj_Jac_Tensor_i[3][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[3][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[3][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[3][2] += factor*val_Proj_Visc_Flux[2];
    
    
	} 
	else {

		su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
		su2double thetaz = theta + val_normal[2]*val_normal[2]/3.0;

		su2double etax = val_normal[1]*val_normal[2]/3.0;
		su2double etay = val_normal[0]*val_normal[2]/3.0;
		su2double etaz = val_normal[0]*val_normal[1]/3.0;

		su2double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz   + val_Mean_PrimVar[3]*etay;
		su2double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay + val_Mean_PrimVar[3]*etax;
		su2double piz = val_Mean_PrimVar[1]*etay   + val_Mean_PrimVar[2]*etax   + val_Mean_PrimVar[3]*thetaz;

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
		val_Proj_Jac_Tensor_i[4][0] = -factor*(rhoovisc*theta*(phi_rho+phi*phi_p) - (pix*val_Mean_PrimVar[1] + piy*val_Mean_PrimVar[2] + piz*val_Mean_PrimVar[3]));
		val_Proj_Jac_Tensor_i[4][1] = -factor*(pix-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[1]);
		val_Proj_Jac_Tensor_i[4][2] = -factor*(piy-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[2]);
		val_Proj_Jac_Tensor_i[4][3] = -factor*(piz-rhoovisc*theta*phi_p*(Gamma-1)*val_Mean_PrimVar[3]);
		val_Proj_Jac_Tensor_i[4][4] = -factor*((Gamma-1)*rhoovisc*theta*phi_p);

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];
    
		factor = 0.5/Density;
		val_Proj_Jac_Tensor_i[4][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[4][0] += factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_i[4][3] += factor*val_Proj_Visc_Flux[3];
		val_Proj_Jac_Tensor_j[4][3] += factor*val_Proj_Visc_Flux[3];



	}

//			for (iVar = 0; iVar < nVar; iVar++) {
//				for (jVar = 0; jVar < nVar; jVar++) {
//					cout << val_Proj_Jac_Tensor_i[iVar][jVar] << " " << val_Proj_Jac_Tensor_j[iVar][jVar] << endl;
//				}
//			}
//	        getchar();

}

void CNumerics::GetViscousProjJacs(su2double *val_Mean_PrimVar,
									su2double **val_gradprimvar,
									su2double *val_Mean_SecVar,
									su2double val_laminar_viscosity,
									su2double val_eddy_viscosity,
									su2double val_thermal_conductivity,
									su2double val_heat_capacity_cp,
									su2double val_dist_ij,
									su2double *val_normal, su2double val_dS,
									su2double *val_Proj_Visc_Flux,
									su2double **val_Proj_Jac_Tensor_i,
									su2double **val_Proj_Jac_Tensor_j) {

  AD_BEGIN_PASSIVE
	/* Viscous flux Jacobians for arbitrary equations of state */

	//order of val_mean_primitives: T, vx, vy, vz, P, rho, ht
	//order of secondary:dTdrho_e, dTde_rho
	unsigned short iDim, iVar, jVar;

	su2double sqvel = 0.0, theta= 0.0, proj_viscousflux_vel= 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		theta += val_normal[iDim]*val_normal[iDim];
		sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
		proj_viscousflux_vel += val_Proj_Visc_Flux[iDim+1]*val_Mean_PrimVar[iDim+1];
	}

	su2double rho = val_Mean_PrimVar[nDim+2];
	su2double P= val_Mean_PrimVar[nDim+1];
	su2double h= val_Mean_PrimVar[nDim+3];
	su2double dTdrho_e= val_Mean_SecVar[0];
	su2double dTde_rho= val_Mean_SecVar[1];


	su2double dTdu0= dTdrho_e + dTde_rho*(-(h-P/rho) + sqvel)*(1/rho);
	su2double dTdu1= dTde_rho*(-val_Mean_PrimVar[1])*(1/rho);
	su2double dTdu2= dTde_rho*(-val_Mean_PrimVar[2])*(1/rho);
	su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	su2double total_conductivity = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;


	su2double factor1 = total_viscosity*val_dS/(rho*val_dist_ij);
	su2double factor2 = total_conductivity*val_dS/val_dist_ij;


	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		for (unsigned short jVar = 0; jVar < nVar; jVar++) {
			val_Proj_Jac_Tensor_i[iVar][jVar] = 0.0;
			val_Proj_Jac_Tensor_j[iVar][jVar] = 0.0;
		}
	}

	if (nDim == 2) {

	    /* 2D Jacobian: (Fv1, Fv2, Fv3, Fv4) --> (T, vx, vy, rho) */

		su2double dTdu3= dTde_rho*(1/rho);

		su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;

		su2double etaz = val_normal[0]*val_normal[1]/3.0;

		su2double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz;
		su2double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = factor1*pix;
		val_Proj_Jac_Tensor_i[1][1] = -factor1*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor1*etaz;
		val_Proj_Jac_Tensor_i[1][3] = 0.0;
		val_Proj_Jac_Tensor_i[2][0] = factor1*piy;
		val_Proj_Jac_Tensor_i[2][1] = -factor1*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor1*thetay;
		val_Proj_Jac_Tensor_i[2][3] = 0.0;

		val_Proj_Jac_Tensor_i[3][0] = val_Proj_Jac_Tensor_i[1][0]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][0]*val_Mean_PrimVar[2];
		val_Proj_Jac_Tensor_i[3][0] += -factor2*theta*dTdu0;
		val_Proj_Jac_Tensor_i[3][1] = val_Proj_Jac_Tensor_i[1][1]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][1]*val_Mean_PrimVar[2];
		val_Proj_Jac_Tensor_i[3][1] += -factor2*theta*dTdu1;
		val_Proj_Jac_Tensor_i[3][2] = val_Proj_Jac_Tensor_i[1][2]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][2]*val_Mean_PrimVar[2];
		val_Proj_Jac_Tensor_i[3][2] += -factor2*theta*dTdu2;
		val_Proj_Jac_Tensor_i[3][3] = -factor2*theta*dTdu3;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		su2double factor = 0.5/rho;
		val_Proj_Jac_Tensor_i[3][0] -= factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[3][0] -= factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[3][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[3][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[3][2] += factor*val_Proj_Visc_Flux[2];




	} else {


		su2double dTdu3= dTde_rho*(-val_Mean_PrimVar[3])*(1/rho);
		su2double dTdu4= dTde_rho*(1/rho);

		su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
		su2double thetaz = theta + val_normal[2]*val_normal[2]/3.0;

		su2double etax = val_normal[1]*val_normal[2]/3.0;
		su2double etay = val_normal[0]*val_normal[2]/3.0;
		su2double etaz = val_normal[0]*val_normal[1]/3.0;

		su2double pix = val_Mean_PrimVar[1]*thetax + val_Mean_PrimVar[2]*etaz   + val_Mean_PrimVar[3]*etay;
		su2double piy = val_Mean_PrimVar[1]*etaz   + val_Mean_PrimVar[2]*thetay + val_Mean_PrimVar[3]*etax;
		su2double piz = val_Mean_PrimVar[1]*etay   + val_Mean_PrimVar[2]*etax   + val_Mean_PrimVar[3]*thetaz;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[0][4] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = factor1*pix;
		val_Proj_Jac_Tensor_i[1][1] = -factor1*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor1*etaz;
		val_Proj_Jac_Tensor_i[1][3] = -factor1*etay;
		val_Proj_Jac_Tensor_i[1][4] = 0.0;
		val_Proj_Jac_Tensor_i[2][0] = factor1*piy;
		val_Proj_Jac_Tensor_i[2][1] = -factor1*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor1*thetay;
		val_Proj_Jac_Tensor_i[2][3] = -factor1*etax;
		val_Proj_Jac_Tensor_i[2][4] = 0.0;
		val_Proj_Jac_Tensor_i[3][0] = factor1*piz;
		val_Proj_Jac_Tensor_i[3][1] = -factor1*etay;
		val_Proj_Jac_Tensor_i[3][2] = -factor1*etax;
		val_Proj_Jac_Tensor_i[3][3] = -factor1*thetaz;
		val_Proj_Jac_Tensor_i[3][4] = 0.0;
		val_Proj_Jac_Tensor_i[4][0] = val_Proj_Jac_Tensor_i[1][0]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][0]*val_Mean_PrimVar[2]+val_Proj_Jac_Tensor_i[3][0]*val_Mean_PrimVar[3];
		val_Proj_Jac_Tensor_i[4][0] +=  -factor2*theta*dTdu0;
		val_Proj_Jac_Tensor_i[4][1] = val_Proj_Jac_Tensor_i[1][1]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][1]*val_Mean_PrimVar[2]+val_Proj_Jac_Tensor_i[3][1]*val_Mean_PrimVar[3];
		val_Proj_Jac_Tensor_i[4][1] +=  -factor2*theta*dTdu1;
		val_Proj_Jac_Tensor_i[4][2] = val_Proj_Jac_Tensor_i[1][2]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][2]*val_Mean_PrimVar[2]+val_Proj_Jac_Tensor_i[3][2]*val_Mean_PrimVar[3];
		val_Proj_Jac_Tensor_i[4][2] +=  -factor2*theta*dTdu2;
		val_Proj_Jac_Tensor_i[4][3] = val_Proj_Jac_Tensor_i[1][3]*val_Mean_PrimVar[1]+val_Proj_Jac_Tensor_i[2][3]*val_Mean_PrimVar[2]+val_Proj_Jac_Tensor_i[3][3]*val_Mean_PrimVar[3];
		val_Proj_Jac_Tensor_i[4][3] +=  -factor2*theta*dTdu3;
		val_Proj_Jac_Tensor_i[4][4] = -factor2*theta*dTdu4;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

		su2double factor = 0.5/rho;
		val_Proj_Jac_Tensor_i[4][0] -= factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_j[4][0] -= factor*proj_viscousflux_vel;
		val_Proj_Jac_Tensor_i[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_j[4][1] += factor*val_Proj_Visc_Flux[1];
		val_Proj_Jac_Tensor_i[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_j[4][2] += factor*val_Proj_Visc_Flux[2];
		val_Proj_Jac_Tensor_i[4][3] += factor*val_Proj_Visc_Flux[3];
		val_Proj_Jac_Tensor_j[4][3] += factor*val_Proj_Visc_Flux[3];




		}


  AD_END_PASSIVE
	}





//void CNumerics::GetViscousProjJacs(su2double *val_Mean_PrimVar,
//									su2double **val_gradprimvar,
//									su2double *val_Mean_SecVar,
//									su2double val_laminar_viscosity,
//									su2double val_eddy_viscosity,
//									su2double val_thermal_conductivity,
//									su2double val_heat_capacity_cp,
//									su2double val_dist_ij,
//									su2double *val_normal, su2double val_dS,
//									su2double *val_Proj_Visc_Flux,
//									su2double **val_Proj_Jac_Tensor_i,
//									su2double **val_Proj_Jac_Tensor_j) {
//
//	/* Viscous flux Jacobians for arbitrary equations of state */
//
//	// order of primitives: T, vx, vy, vz, P, rho, h, c, MuLam, MuEddy, kt, Cp
//	// order of secondary: dPdrho_e, dPde_rho, dTdrho_e, dTde_rho, dmudrho_T, dmudT_rho, dktdrho_T, dktdT_rho
//
//	unsigned short iDim, iVar, jVar;
//	su2double **val_Proj_Jac_Tensor_i_P, **val_Proj_Jac_Tensor_j_P;
//	su2double **val_Jac_PC;
//
//	su2double sqvel = 0.0;
//	for (iDim = 0; iDim < nDim; iDim++) {
//		sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
//	}
//
//	su2double vx = val_Mean_PrimVar[1];
//	su2double vy = val_Mean_PrimVar[2];
//	su2double vz = val_Mean_PrimVar[3];
//	su2double rho = val_Mean_PrimVar[nDim+2];
//	su2double P = val_Mean_PrimVar[nDim+1];
//	su2double dmudrho_T = val_Mean_SecVar[4];
//	su2double dmudT_rho = val_Mean_SecVar[5];
//	su2double dktdrho_T = val_Mean_SecVar[6];
//	su2double dktdT_rho = val_Mean_SecVar[7];
//
//	su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
//	su2double total_conductivity = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;
//
//	val_Proj_Jac_Tensor_i_P = new su2double* [nVar];
//	val_Proj_Jac_Tensor_j_P = new su2double* [nVar];
//	val_Jac_PC = new su2double* [nVar];
//
//	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
//		val_Proj_Jac_Tensor_i_P[iVar] = new su2double [nVar];
//		val_Proj_Jac_Tensor_j_P[iVar] = new su2double [nVar];
//		val_Jac_PC[iVar] = new su2double [nVar];
//	}
//
//	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
//		for (unsigned short jVar = 0; jVar < nVar; jVar++) {
//			val_Proj_Jac_Tensor_i[iVar][jVar] = 0.0;
//			val_Proj_Jac_Tensor_j[iVar][jVar] = 0.0;
//			val_Proj_Jac_Tensor_i_P[iVar][jVar] = 0.0;
//			val_Proj_Jac_Tensor_j_P[iVar][jVar] = 0.0;
//			val_Jac_PC[iVar][jVar] = 0.0;
//		}
//	}
//
//	if (nDim == 2) {
//
//	    /* 2D Jacobian: (Fv1, Fv2, Fv3, Fv4) --> (T, vx, vy, rho) */
//
//		su2double factor1 = 4.0/3.0*pow(val_normal[0],2) + pow(val_normal[1],2);
//		su2double factor2 = 1.0/3.0*val_normal[0]*val_normal[1];
//		su2double factor3 = 4.0/3.0*pow(val_normal[1],2) + pow(val_normal[0],2);
//
//		val_Proj_Jac_Tensor_i_P[0][0] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][1] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][2] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][3] = 0.0;
//
//		val_Proj_Jac_Tensor_i_P[1][0] = 0.5*dmudT_rho*val_Proj_Visc_Flux[1]/total_viscosity;
//		val_Proj_Jac_Tensor_i_P[1][1] = -total_viscosity/val_dist_ij*factor1;
//		val_Proj_Jac_Tensor_i_P[1][2] = -total_viscosity/val_dist_ij*factor2;
//		val_Proj_Jac_Tensor_i_P[1][3] = 0.5*dmudrho_T*val_Proj_Visc_Flux[1]/total_viscosity;
//
//		val_Proj_Jac_Tensor_i_P[2][0] = 0.5*dmudT_rho*val_Proj_Visc_Flux[2]/total_viscosity;
//		val_Proj_Jac_Tensor_i_P[2][1] = -total_viscosity/val_dist_ij*factor2;
//		val_Proj_Jac_Tensor_i_P[2][2] = -total_viscosity/val_dist_ij*factor3;
//		val_Proj_Jac_Tensor_i_P[2][3] = 0.5*dmudrho_T*val_Proj_Visc_Flux[2]/total_viscosity;
//
//		val_Proj_Jac_Tensor_i_P[3][0] = vx*val_Proj_Jac_Tensor_i_P[1][0] + vy*val_Proj_Jac_Tensor_i_P[2][0];
//		val_Proj_Jac_Tensor_i_P[3][1] = vx*val_Proj_Jac_Tensor_i_P[1][1] + vy*val_Proj_Jac_Tensor_i_P[2][1];
//		val_Proj_Jac_Tensor_i_P[3][2] = vx*val_Proj_Jac_Tensor_i_P[1][2] + vy*val_Proj_Jac_Tensor_i_P[2][2];
//		val_Proj_Jac_Tensor_i_P[3][3] = vx*val_Proj_Jac_Tensor_i_P[1][3] + vy*val_Proj_Jac_Tensor_i_P[2][3];
//
//		su2double etax = pow(val_normal[0],2)/val_dist_ij;
//		su2double etay = pow(val_normal[1],2)/val_dist_ij;
//		val_Proj_Jac_Tensor_i_P[3][0] += -total_conductivity*etax + 0.5*val_normal[0]*dktdT_rho*val_gradprimvar[0][0] - total_conductivity*etay + 0.5*val_normal[1]*dktdT_rho*val_gradprimvar[0][1];
//		val_Proj_Jac_Tensor_i_P[3][1] += 1.0/2.0*val_Proj_Visc_Flux[1];
//		val_Proj_Jac_Tensor_i_P[3][2] += 1.0/2.0*val_Proj_Visc_Flux[2];
//		val_Proj_Jac_Tensor_i_P[3][3] += 0.5*val_normal[0]*dktdrho_T*val_gradprimvar[0][0] +0.5*val_normal[1]*dktdrho_T*val_gradprimvar[0][1];
//
//		for (iVar = 0; iVar < nVar; iVar++) {
//			for (jVar = 0; jVar < nVar; jVar++) {
//				val_Proj_Jac_Tensor_i_P[iVar][jVar] *= val_dS;
//		        val_Proj_Jac_Tensor_j_P[iVar][jVar]  = -val_Proj_Jac_Tensor_i_P[iVar][jVar];  /*Jacobian j*/
//			}
//		}
//
//	    /* 2D Jacobian: (T, vx, vy, rho) --> (u1, u2, u3, u4) */
//		GetPrimitive2Conservative (val_Mean_PrimVar, val_Mean_SecVar, val_Jac_PC);
//
//	    /* 2D Jacobian: (Fv1, Fv2, Fv3, Fv4) --> (u1, u2, u3, u4) */
//		for (iVar = 0; iVar < nVar; iVar++) {
//			for (jVar = 0; jVar < nVar; jVar++) {
//				for (unsigned short kVar = 0; kVar < nVar; kVar++) {
//					val_Proj_Jac_Tensor_i[iVar][jVar] += val_Proj_Jac_Tensor_i_P[iVar][kVar]*val_Jac_PC[kVar][jVar];
//					val_Proj_Jac_Tensor_j[iVar][jVar] += val_Proj_Jac_Tensor_j_P[iVar][kVar]*val_Jac_PC[kVar][jVar];
//			    }
//			}
//		}
//
////		for (iVar = 0; iVar < nVar; iVar++) {
////			for (jVar = 0; jVar < nVar; jVar++) {
////				cout << val_Proj_Jac_Tensor_i[iVar][jVar] << " " << val_Proj_Jac_Tensor_j[iVar][jVar] << endl;
////			}
////		}
////        getchar();
//
//	}
//	else {
//
//	    /* 3D Jacobian: (Fv1, Fv2, Fv3, Fv4, Fv5) --> (T, vx, vy, vz, rho) */
//
////		su2double factor1 = 4.0/3.0*pow(val_normal[0],2) + pow(val_normal[1],2) + pow(val_normal[2],2);
////		su2double factor2 = 1.0/3.0*val_normal[0]*val_normal[1];
////		su2double factor3 = 1.0/3.0*val_normal[0]*val_normal[2];
////		su2double factor4 = pow(val_normal[0],2) + 4.0/3.0*pow(val_normal[1],2) + pow(val_normal[2],2);
////        su2double factor5 = 1.0/3.0*val_normal[1]*val_normal[2];
////        su2double factor6 = pow(val_normal[0],2) + pow(val_normal[1],2) + 4.0/3.0*pow(val_normal[2],2);
//
//		su2double thetax = 4.0/3.0*pow(val_normal[0],2) + pow(val_normal[1],2) + pow(val_normal[2],2);
//		su2double thetay = pow(val_normal[0],2) + 4.0/3.0*pow(val_normal[1],2) + pow(val_normal[2],2);
//        su2double thetaz = pow(val_normal[0],2) + pow(val_normal[1],2) + 4.0/3.0*pow(val_normal[2],2);
//		su2double pix = 1.0/3.0*val_normal[0]*val_normal[1];
//		su2double piy = 1.0/3.0*val_normal[0]*val_normal[2];
//		su2double piz = 1.0/3.0*val_normal[1]*val_normal[2];
//
//		val_Proj_Jac_Tensor_i_P[0][0] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][1] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][2] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][3] = 0.0;
//		val_Proj_Jac_Tensor_i_P[0][4] = 0.0;
//
//		val_Proj_Jac_Tensor_i_P[1][0] = 0.5*dmudT_rho*val_Proj_Visc_Flux[1]/total_viscosity;
//		val_Proj_Jac_Tensor_i_P[1][1] = -total_viscosity/val_dist_ij*thetax;
//		val_Proj_Jac_Tensor_i_P[1][2] = -total_viscosity/val_dist_ij*pix;
//		val_Proj_Jac_Tensor_i_P[1][3] = -total_viscosity/val_dist_ij*piy;
//		val_Proj_Jac_Tensor_i_P[1][4] = 0.5*dmudrho_T*val_Proj_Visc_Flux[1]/total_viscosity;
//
//		val_Proj_Jac_Tensor_i_P[2][0] = 0.5*dmudT_rho*val_Proj_Visc_Flux[2]/total_viscosity;
//		val_Proj_Jac_Tensor_i_P[2][1] = -total_viscosity/val_dist_ij*pix;
//		val_Proj_Jac_Tensor_i_P[2][2] = -total_viscosity/val_dist_ij*thetay;
//		val_Proj_Jac_Tensor_i_P[2][3] = -total_viscosity/val_dist_ij*piz;
//		val_Proj_Jac_Tensor_i_P[2][4] = 0.5*dmudrho_T*val_Proj_Visc_Flux[2]/total_viscosity;
//
//		val_Proj_Jac_Tensor_i_P[3][0] = 0.5*dmudT_rho*val_Proj_Visc_Flux[3]/total_viscosity;
//		val_Proj_Jac_Tensor_i_P[3][1] = -total_viscosity/val_dist_ij*piy;
//		val_Proj_Jac_Tensor_i_P[3][2] = -total_viscosity/val_dist_ij*piz;
//		val_Proj_Jac_Tensor_i_P[3][3] = -total_viscosity/val_dist_ij*thetaz;
//		val_Proj_Jac_Tensor_i_P[3][4] = 0.5*dmudrho_T*val_Proj_Visc_Flux[3]/total_viscosity;
//
//		val_Proj_Jac_Tensor_i_P[4][0] = vx*val_Proj_Jac_Tensor_i_P[1][0] + vy*val_Proj_Jac_Tensor_i_P[2][0] + vz*val_Proj_Jac_Tensor_i_P[3][0];
//		val_Proj_Jac_Tensor_i_P[4][1] = vx*val_Proj_Jac_Tensor_i_P[1][1] + vy*val_Proj_Jac_Tensor_i_P[2][1] + vz*val_Proj_Jac_Tensor_i_P[3][1];
//		val_Proj_Jac_Tensor_i_P[4][2] = vx*val_Proj_Jac_Tensor_i_P[1][2] + vy*val_Proj_Jac_Tensor_i_P[2][2] + vz*val_Proj_Jac_Tensor_i_P[3][2];
//		val_Proj_Jac_Tensor_i_P[4][3] = vx*val_Proj_Jac_Tensor_i_P[1][3] + vy*val_Proj_Jac_Tensor_i_P[2][3] + vz*val_Proj_Jac_Tensor_i_P[3][3];
//		val_Proj_Jac_Tensor_i_P[4][4] = vx*val_Proj_Jac_Tensor_i_P[1][4] + vy*val_Proj_Jac_Tensor_i_P[2][4] + vz*val_Proj_Jac_Tensor_i_P[3][4];
//
//		su2double etax = pow(val_normal[0],2)/val_dist_ij;
//		su2double etay = pow(val_normal[1],2)/val_dist_ij;
//		su2double etaz = pow(val_normal[2],2)/val_dist_ij;
//		val_Proj_Jac_Tensor_i_P[4][0] += -total_conductivity*etax - total_conductivity*etay - total_conductivity*etaz;
//		val_Proj_Jac_Tensor_i_P[4][1] += 1.0/2.0*val_Proj_Visc_Flux[1];
//		val_Proj_Jac_Tensor_i_P[4][2] += 1.0/2.0*val_Proj_Visc_Flux[2];
//		val_Proj_Jac_Tensor_i_P[4][3] += 1.0/2.0*val_Proj_Visc_Flux[3];
//
//		val_Proj_Jac_Tensor_i_P[4][0] += 0.5*val_normal[0]*dktdT_rho*val_gradprimvar[0][0] + 0.5*val_normal[1]*dktdT_rho*val_gradprimvar[0][1] + 0.5*val_normal[2]*dktdT_rho*val_gradprimvar[0][2];
//		val_Proj_Jac_Tensor_i_P[4][4] += 0.5*val_normal[0]*dktdrho_T*val_gradprimvar[0][0] + 0.5*val_normal[1]*dktdrho_T*val_gradprimvar[0][1] + 0.5*val_normal[2]*dktdrho_T*val_gradprimvar[0][2];
//
//		for (iVar = 0; iVar < nVar; iVar++) {
//			for (jVar = 0; jVar < nVar; jVar++) {
//				val_Proj_Jac_Tensor_i_P[iVar][jVar] *= val_dS;
//				val_Proj_Jac_Tensor_j_P[iVar][jVar] = -val_Proj_Jac_Tensor_i_P[iVar][jVar];	 /*Jacobian j*/
//			}
//		}
//
//	    /* 3D Jacobian: (T, vx, vy, vz, rho) --> (u1, u2, u3, u4, u5) */
//		GetPrimitive2Conservative (val_Mean_PrimVar, val_Mean_SecVar, val_Jac_PC);
//
//	    /* 3D Jacobian: (Fv1, Fv2, Fv3, Fv4, Fv5) --> (u1, u2, u3, u4, u5) */
//		for (iVar = 0; iVar < nVar; iVar++) {
//			for (jVar = 0; jVar < nVar; jVar++) {
//				for (unsigned short kVar = 0; kVar < nVar; kVar++) {
//					val_Proj_Jac_Tensor_i[iVar][jVar] += val_Proj_Jac_Tensor_i_P[iVar][kVar]*val_Jac_PC[kVar][jVar];
//					val_Proj_Jac_Tensor_j[iVar][jVar] += val_Proj_Jac_Tensor_j_P[iVar][kVar]*val_Jac_PC[kVar][jVar];
//			    }
//			}
//		}
//
////				for (iVar = 0; iVar < nVar; iVar++) {
////					for (jVar = 0; jVar < nVar; jVar++) {
////						cout << val_Proj_Jac_Tensor_i[iVar][jVar] << " " << val_Proj_Jac_Tensor_j[iVar][jVar] << endl;
////					}
////				}
////		        getchar();
//
//	}
//
//    /*--- Deallocate ---*/
//	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
//		delete [] val_Proj_Jac_Tensor_i_P[iVar];
//		delete [] val_Proj_Jac_Tensor_j_P[iVar];
//		delete [] val_Jac_PC[iVar];
//	}
//	delete [] val_Proj_Jac_Tensor_i_P;
//	delete [] val_Proj_Jac_Tensor_j_P;
//	delete [] val_Jac_PC;
//
//}

void CNumerics::GetPrimitive2Conservative (su2double *val_Mean_PrimVar, su2double *val_Mean_SecVar, su2double **val_Jac_PC) {

	unsigned short iVar, jVar, iDim;

	// order of primitives: T, vx, vy, vz, P, rho, h, c, MuLam, MuEddy, kt, Cp
	// order of secondary: dPdrho_e, dPde_rho, dTdrho_e, dTde_rho, dmudrho_T, dmudT_rho, dktdrho_T, dktdT_rho

	su2double vx = val_Mean_PrimVar[1];
	su2double vy = val_Mean_PrimVar[2];
	su2double vz = val_Mean_PrimVar[3];
	su2double rho = val_Mean_PrimVar[nDim+2];
	su2double P = val_Mean_PrimVar[nDim+1];
	su2double e = val_Mean_PrimVar[nDim+3] - P/rho;
	su2double dTdrho_e = val_Mean_SecVar[2];
	su2double dTde_rho = val_Mean_SecVar[3];

	su2double sqvel = 0.0;

	for (iDim = 0; iDim < nDim; iDim++) {
		sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
	}

  /*--- Initialize the Jacobian matrix ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jac_PC[iVar][jVar] = 0.0;
    }
  }

  /*--- Primitives to conservatives Jacobian matrix : (T, vx, vy, vz, rho) --> (u1, u2, u3, u4, u5) ---*/
  if (nDim == 2) {

	val_Jac_PC[0][0] = dTdrho_e - e/rho*dTde_rho + 0.5*dTde_rho*sqvel/rho;
	val_Jac_PC[0][1] = -1/rho*dTde_rho*vx;
	val_Jac_PC[0][2] = -1/rho*dTde_rho*vy;
	val_Jac_PC[0][3] = 1/rho*dTde_rho;

	val_Jac_PC[1][0] = -vx/rho;
	val_Jac_PC[1][1] = 1/rho;
	val_Jac_PC[1][2] = 0.0;
	val_Jac_PC[1][3] = 0.0;

	val_Jac_PC[2][0] = -vy/rho;
	val_Jac_PC[2][1] = 0.0;
	val_Jac_PC[2][2] = 1/rho;
	val_Jac_PC[2][3] = 0.0;

	val_Jac_PC[3][0] = 1.0;
	val_Jac_PC[3][1] = 0.0;
	val_Jac_PC[3][2] = 0.0;
	val_Jac_PC[3][3] = 0.0;

  }
  else {

	val_Jac_PC[0][0] = dTdrho_e - e/rho*dTde_rho + 0.5*dTde_rho*sqvel/rho;
	val_Jac_PC[0][1] = -1/rho*dTde_rho*vx;
	val_Jac_PC[0][2] = -1/rho*dTde_rho*vy;
	val_Jac_PC[0][3] = -1/rho*dTde_rho*vz;
	val_Jac_PC[0][4] = 1/rho*dTde_rho;

	val_Jac_PC[1][0] = -vx/rho;
	val_Jac_PC[1][1] = 1/rho;
	val_Jac_PC[1][2] = 0.0;
	val_Jac_PC[1][3] = 0.0;
	val_Jac_PC[1][4] = 0.0;

	val_Jac_PC[2][0] = -vy/rho;
	val_Jac_PC[2][1] = 0.0;
	val_Jac_PC[2][2] = 1/rho;
	val_Jac_PC[2][3] = 0.0;
	val_Jac_PC[2][4] = 0.0;

	val_Jac_PC[3][0] = -vz/rho;
	val_Jac_PC[3][1] = 0.0;
	val_Jac_PC[3][2] = 0.0;
	val_Jac_PC[3][3] = 1/rho;
	val_Jac_PC[3][4] = 0.0;

	val_Jac_PC[4][0] = 1.0;
	val_Jac_PC[4][1] = 0.0;
	val_Jac_PC[4][2] = 0.0;
	val_Jac_PC[4][3] = 0.0;
	val_Jac_PC[4][4] = 0.0;

  }
}

void CNumerics::GetViscousArtCompProjJacs(su2double val_laminar_viscosity,
		su2double val_eddy_viscosity, su2double val_dist_ij, su2double *val_normal, su2double val_dS,
		su2double **val_Proj_Jac_Tensor_i, su2double **val_Proj_Jac_Tensor_j) {
	unsigned short iDim, iVar, jVar;

	su2double theta = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		theta += val_normal[iDim]*val_normal[iDim];

	su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
	su2double factor = total_viscosity/(val_dist_ij)*val_dS;

	if (nDim == 3) {
		su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
		su2double thetaz = theta + val_normal[2]*val_normal[2]/3.0;

		su2double etax = val_normal[1]*val_normal[2]/3.0;
		su2double etay = val_normal[0]*val_normal[2]/3.0;
		su2double etaz = val_normal[0]*val_normal[1]/3.0;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[0][3] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = 0.0;
		val_Proj_Jac_Tensor_i[1][1] = -factor*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor*etaz;
		val_Proj_Jac_Tensor_i[1][3] = -factor*etay;
		val_Proj_Jac_Tensor_i[2][0] = 0.0;
		val_Proj_Jac_Tensor_i[2][1] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor*thetay;
		val_Proj_Jac_Tensor_i[2][3] = -factor*etax;
		val_Proj_Jac_Tensor_i[3][0] = 0.0;
		val_Proj_Jac_Tensor_i[3][1] = -factor*etay;
		val_Proj_Jac_Tensor_i[3][2] = -factor*etax;
		val_Proj_Jac_Tensor_i[3][3] = -factor*thetaz;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

	}

	if (nDim == 2) {
		su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
		su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
		su2double etaz = val_normal[0]*val_normal[1]/3.0;

		val_Proj_Jac_Tensor_i[0][0] = 0.0;
		val_Proj_Jac_Tensor_i[0][1] = 0.0;
		val_Proj_Jac_Tensor_i[0][2] = 0.0;
		val_Proj_Jac_Tensor_i[1][0] = 0.0;
		val_Proj_Jac_Tensor_i[1][1] = -factor*thetax;
		val_Proj_Jac_Tensor_i[1][2] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][0] = 0.0;
		val_Proj_Jac_Tensor_i[2][1] = -factor*etaz;
		val_Proj_Jac_Tensor_i[2][2] = -factor*thetay;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];
	}
  
}

void CNumerics::CreateBasis(su2double *val_Normal) {
  unsigned short iDim;
  su2double modm, modl;
  
  /*--- Define l as a vector in the plane normal to the supplied vector ---*/
  l[0] = 0.0;
  l[1] = -val_Normal[2];
  l[2] = val_Normal[1];
  
  /*--- Check for the zero vector and re-assign if needed ---*/
  if (l[0] == 0.0 && l[1] == 0.0 && l[2] == 0.0) {
    l[0] = -val_Normal[2];
    l[1] = 0.0;
    l[2] = val_Normal[0];
  }
  
  /*--- Take vector product of n * l to make m ---*/
  m[0] = val_Normal[1]*l[2] - val_Normal[2]*l[1];
	m[1] = val_Normal[2]*l[0] - val_Normal[0]*l[2];
	m[2] = val_Normal[0]*l[1] - val_Normal[1]*l[0];
  
  /*--- Normalize ---*/
  modm =0 ; modl = 0;
  for (iDim =0 ; iDim < nDim; iDim++) {
    modm += m[iDim]*m[iDim];
    modl += l[iDim]*l[iDim];
  }
  modm = sqrt(modm);
  modl = sqrt(modl);
  for (iDim =0 ; iDim < nDim; iDim++) {
    l[iDim] = l[iDim]/modl;
    m[iDim] = m[iDim]/modm;
  }
}

CSourceNothing::CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceNothing::~CSourceNothing(void) { }
