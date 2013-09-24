/*!
 * \file numerics_adjoint_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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
#include <limits>

CUpwRoe_AdjTNE2::CUpwRoe_AdjTNE2(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim,
                                                              val_nVar,
                                                              config) {
  
  unsigned short iVar;
  
  /*--- Read configuration parameters ---*/
	implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  ionization = (config->GetIonization());
  
  /*--- Define useful constants ---*/
  nVar = val_nVar;
  nDim = val_nDim;
  nSpecies = config->GetnSpecies();
  
  /*--- Allocate arrays ---*/
	Residual_Roe       = new double [nVar];
  l                  = new double [nDim];
  m                  = new double [nDim];
	RoeVelocity        = new double [nDim];
	Velocity_i         = new double [nDim];
	Velocity_j         = new double [nDim];
  RoeDensity         = new double [nSpecies];
  Density_i          = new double [nSpecies];
  Density_j          = new double [nSpecies];
  dPdrhos            = new double [nSpecies];
	Lambda             = new double [nVar];
	P_Tensor           = new double* [nVar];
	invP_Tensor        = new double* [nVar];
	Proj_flux_tensor_i = new double* [nVar];
	Proj_flux_tensor_j = new double* [nVar];
  Proj_Jac_Tensor_i  = new double* [nVar];
	Proj_Jac_Tensor_j  = new double* [nVar];
	Proj_ModJac_Tensor = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]           = new double [nVar];
		invP_Tensor[iVar]        = new double [nVar];
    Proj_Jac_Tensor_i[iVar]  = new double [nVar];
    Proj_Jac_Tensor_j[iVar]  = new double [nVar];
		Proj_flux_tensor_i[iVar] = new double [nVar];
		Proj_flux_tensor_j[iVar] = new double [nVar];
		Proj_ModJac_Tensor[iVar] = new double [nVar];
	}
  
}

CUpwRoe_AdjTNE2::~CUpwRoe_AdjTNE2(void) {
  unsigned short iVar;
  
	delete [] Residual_Roe;
	delete [] RoeVelocity;
	delete [] Velocity_i;
	delete [] Velocity_j;
  delete [] RoeDensity;
  delete [] Density_i;
  delete [] Density_j;
  delete [] dPdrhos;
	delete [] Lambda;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] Proj_flux_tensor_i[iVar];
		delete [] Proj_flux_tensor_j[iVar];
    delete [] Proj_Jac_Tensor_i[iVar];
    delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
  delete [] Proj_Jac_Tensor_i;
  delete [] Proj_Jac_Tensor_j;
	delete [] Proj_ModJac_Tensor;
  delete [] l;
  delete [] m;
  
}

void CUpwRoe_AdjTNE2::ComputeResidual (double *val_residual_i,
                                       double *val_residual_j,
                                       double **val_Jacobian_ii,
                                       double **val_Jacobian_ij,
                                       double **val_Jacobian_ji,
                                       double **val_Jacobian_jj,
                                       CConfig *config) {
  
  unsigned short iDim, iSpecies, iVar, jVar, kVar, nHeavy, nEl;
  double dPdrhoE_i, dPdrhoE_j, dPdrhoEve_i, dPdrhoEve_j;
  double rho_el, rho_el_i, rho_el_j, conc, conc_i, conc_j;
  double Ru, *Ms, rhoCvtr_i, rhoCvtr_j, rhoCvve_i, rhoCvve_j;
  double ProjVelocity_i, ProjVelocity_j, ProjVelocity;
  double dPdrhoE, dPdrhoEve;
  
	/*--- Compute face area ---*/
	Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
  /*--- Compute and unit normal vector ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Determine the number of heavy particle species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {    
    Density_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    Density_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = V_i[VEL_INDEX+iDim];
		Velocity_j[iDim] = V_j[VEL_INDEX+iDim];
	}
	Pressure_i   = V_i[P_INDEX];    Pressure_j   = V_j[P_INDEX];
  Enthalpy_i   = V_i[H_INDEX];    Enthalpy_j   = V_j[H_INDEX];
  SoundSpeed_i = V_i[A_INDEX];    SoundSpeed_j = V_j[A_INDEX];
  Energy_ve_i  = U_i[nSpecies+nDim+1] / V_i[RHO_INDEX];
  Energy_ve_j  = U_j[nSpecies+nDim+1] / V_j[RHO_INDEX];
  rhoCvtr_i    = V_i[RHOCVTR_INDEX];
  rhoCvtr_j    = V_j[RHOCVTR_INDEX];
  rhoCvve_i    = V_i[RHOCVVE_INDEX];
  rhoCvve_j    = V_j[RHOCVVE_INDEX];
  
  /*--- Calculate mean quantities ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    RoeDensity[iSpecies] = 0.5 * (Density_i[iSpecies] + Density_j[iSpecies]);
  for (iDim = 0; iDim < nDim; iDim++)
    RoeVelocity[iDim] = 0.5 * (Velocity_i[iDim] + Velocity_j[iDim]);
  RoeEnthalpy   = 0.5 * (Enthalpy_i   + Enthalpy_j);
  RoeEnergy_ve  = 0.5 * (Energy_ve_i  + Energy_ve_j);
  RoeSoundSpeed = 0.5 * (SoundSpeed_i + SoundSpeed_j);
  
  /*--- Read parameters from CConfig ---*/
  Ms = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  
  /*--- Calculate useful quantities ---*/
  if (ionization) {
    rho_el   = RoeDensity[nSpecies-1];
    rho_el_i = Density_i[nSpecies-1];
    rho_el_j = Density_j[nSpecies-1];
  } else {
    rho_el   = 0.0;
    rho_el_i = 0.0;
    rho_el_j = 0.0;
  }
  conc   = 0.0;
  conc_i = 0.0;
  conc_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc   += RoeDensity[iSpecies] / Ms[iSpecies];
    conc_i += Density_i[iSpecies] / Ms[iSpecies];
    conc_j += Density_j[iSpecies] / Ms[iSpecies];
    dPdrhos[iSpecies] = 0.5*(dPdrhos_i[iSpecies] + dPdrhos_j[iSpecies]);
  }
  dPdrhoE_i   = conc_i*Ru / rhoCvtr_i;
  dPdrhoE_j   = conc_j*Ru / rhoCvtr_j;
  dPdrhoEve_i = -dPdrhoE_i + rho_el_i * Ru/Ms[nSpecies-1] * 1.0/rhoCvve_i;
  dPdrhoEve_j = -dPdrhoE_j + rho_el_j * Ru/Ms[nSpecies-1] * 1.0/rhoCvve_j;
  dPdrhoE     = conc*Ru / (0.5*(rhoCvtr_i+rhoCvtr_j));
  dPdrhoEve   = -dPdrhoE + rho_el*Ru/Ms[nSpecies-1]
                           *0.5*(rhoCvve_i+rhoCvve_j);

	/*--- Calculate Projected Flux Jacobians (inviscid) ---*/
  // Note: Scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal
  GetInviscidProjJac(Density_i, Velocity_i, &Enthalpy_i, &Energy_ve_i,
                     dPdrhos_i, dPdrhoE_i, dPdrhoEve_i, Normal, 0.5,
                     Proj_Jac_Tensor_i);
  GetInviscidProjJac(Density_j, Velocity_j, &Enthalpy_j, &Energy_ve_j,
                     dPdrhos_j, dPdrhoE_j, dPdrhoEve_j, Normal, 0.5,
                     Proj_Jac_Tensor_j);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar] = 0.0; val_residual_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
			val_residual_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
		}
	}
  
  /*--- Calculate dual grid tangent vectors for P & invP ---*/
  CreateBasis(UnitNormal);
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix(RoeDensity, RoeVelocity, &RoeEnthalpy, &RoeEnergy_ve,
             &RoeSoundSpeed, dPdrhos, dPdrhoE, dPdrhoEve, UnitNormal,
             l, m, P_Tensor);
  GetPMatrix_inv(RoeDensity, RoeVelocity, &RoeEnergy_ve, &RoeSoundSpeed,
                 dPdrhos, dPdrhoE, dPdrhoEve, UnitNormal, l, m, invP_Tensor);
  
  /*--- Compute projected velocities ---*/
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Calculate eigenvalues ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Lambda[iSpecies] = ProjVelocity;
  for (iDim = 0; iDim < nDim-1; iDim++)
    Lambda[nSpecies+iDim] = ProjVelocity;
  Lambda[nSpecies+nDim-1] = ProjVelocity + RoeSoundSpeed;
  Lambda[nSpecies+nDim]   = ProjVelocity - RoeSoundSpeed;
  Lambda[nSpecies+nDim+1] = ProjVelocity;
  
	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);
  
	/*--- Flux approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = 0.5*Proj_ModJac_Tensor_ij;
		}
	}
  
	for (iVar = 0; iVar < nVar; iVar++)
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] -= Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
			val_residual_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
		}
  
	/*--- Implicit contributions, Transpose the matrices and store the Jacobians. Note the negative
	 sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_ii[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ij[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ji[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
				val_Jacobian_jj[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
			}
		}
	}
}


void CUpwRoe_AdjTNE2::CreateBasis(double *val_Normal) {
  unsigned short iDim;
  double modm, modl;
  
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




CCentJST_AdjTNE2::CCentJST_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
  
	Diff_Psi = new double [nVar]; Diff_Lapl = new double [nVar];
	Und_Lapl_i = new double [nVar]; Und_Lapl_j = new double [nVar];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
	MeanPhi = new double [nDim];
  
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_AdjTNE2();
	Param_Kappa_4 = config->GetKappa_4th_AdjTNE2();
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
}

CCentJST_AdjTNE2::~CCentJST_AdjTNE2(void) {
  
	delete [] Diff_Psi; delete [] Diff_Lapl;
	delete [] Und_Lapl_i; delete [] Und_Lapl_j;
	delete [] Velocity_i; delete [] Velocity_j;
	delete [] MeanPhi;
}

void CCentJST_AdjTNE2::ComputeResidual (double *val_resconv_i, double *val_resvisc_i,
                                        double *val_resconv_j, double *val_resvisc_j,
                                        double **val_Jacobian_ii, double **val_Jacobian_ij,
                                        double **val_Jacobian_ji, double **val_Jacobian_jj,
                                        CConfig *config) {
  
	/*--- Mean Values ---*/
	MeanPsiRho =  0.5*(Psi_i[0]+Psi_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanPhi[iDim] =  0.5*(Psi_i[iDim+1]+Psi_j[iDim+1]);
	MeanPsiE =  0.5*(Psi_i[nVar-1]+Psi_j[nVar-1]);
  
	/*--- Point i convective residual evaluation ---*/
	ProjVelocity_i = 0; ProjPhi = 0; ProjPhi_Vel = 0; sq_vel = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / U_i[0];
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjPhi += MeanPhi[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_i[iDim];
		sq_vel += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	phis1 = ProjPhi + ProjVelocity_i*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_i*MeanPsiE;
  
	val_resconv_i[0] = ProjVelocity_i*MeanPsiRho - phis2*ProjVelocity_i + Gamma_Minus_One*phis1*sq_vel;
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_i[iDim+1] = ProjVelocity_i*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_i[iDim];
	val_resconv_i[nVar-1] = ProjVelocity_i*MeanPsiE + Gamma_Minus_One*phis1;
  
	/*--- Flux contributions due to grid movement at point i (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_i[0] -= ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjGridVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjGridVel*MeanPsiE;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[0][jDim+1] = -0.5*ProjVelocity_i*Velocity_i[jDim] + Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_ii[0][nVar-1] = 0.5*ProjVelocity_i*(Gamma_Minus_One*sq_vel - Enthalpy_i);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_ii[iDim+1][0] = 0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_ii[iDim+1][jDim+1] = 0.5*Normal[iDim]*Velocity_i[jDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*Normal[jDim];
			val_Jacobian_ii[iDim+1][iDim+1] += 0.5*ProjVelocity_i;
			val_Jacobian_ii[iDim+1][nVar-1] = 0.5*Enthalpy_i*Normal[iDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*ProjVelocity_i;
		}
		val_Jacobian_ii[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[nVar-1][jDim+1] = 0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_ii[nVar-1][nVar-1] = 0.5*Gamma*ProjVelocity_i;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = val_Jacobian_ii[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point i (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
			}
		}
	}
  
  
	/*--- Point j convective residual evaluation ---*/
	ProjVelocity_j = 0; ProjPhi_Vel = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / U_j[0];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_j[iDim];
		sq_vel += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
  
	phis1 = ProjPhi + ProjVelocity_j*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_j*MeanPsiE;
  
	val_resconv_j[0] = -(ProjVelocity_j*MeanPsiRho - phis2*ProjVelocity_j + Gamma_Minus_One*phis1*sq_vel);
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_j[iDim+1] = -(ProjVelocity_j*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_j[iDim]);
	val_resconv_j[nVar-1] = -(ProjVelocity_j*MeanPsiE + Gamma_Minus_One*phis1);
  
	/*--- Flux contributions due to grid movement at point j (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_j[0] += ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjGridVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjGridVel*MeanPsiE;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_jj[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[0][jDim+1] = 0.5*ProjVelocity_j*Velocity_j[jDim] - Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_jj[0][nVar-1] = -0.5*ProjVelocity_j*(Gamma_Minus_One*sq_vel - Enthalpy_j);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_jj[iDim+1][0] = -0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_jj[iDim+1][jDim+1] = -0.5*Normal[iDim]*Velocity_j[jDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*Normal[jDim];
			val_Jacobian_jj[iDim+1][iDim+1] -= 0.5*ProjVelocity_j;
			val_Jacobian_jj[iDim+1][nVar-1] = -0.5*Enthalpy_j*Normal[iDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*ProjVelocity_j;
		}
		val_Jacobian_jj[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[nVar-1][jDim+1] = -0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_jj[nVar-1][nVar-1] = -0.5*Gamma*ProjVelocity_j;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = val_Jacobian_jj[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point j (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
	}
  
	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
	}
  
	/*--- Adjustment to projected velocity due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j += ProjGridVel;
	}
  
	/*--- Compute the spectral radius and stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	double sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
	/*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
		if (implicit) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
		}
	}
}


CCentLax_AdjTNE2::CCentLax_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Normal_ij  = new double [nDim];
  Normal_ji  = new double [nDim];
	DiffPsi   = new double [nVar];
  MeanPsi    = new double [nVar];
  MeanPsiRho = new double [nSpecies];
  MeanPhi    = new double [nDim];
	Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  Density_i = new double [nSpecies];
  Density_j = new double [nSpecies];
  
  Proj_Jac_Tensor_i = new double*[nVar];
  Proj_Jac_Tensor_j = new double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Proj_Jac_Tensor_i[iVar] = new double [nVar];
    Proj_Jac_Tensor_j[iVar] = new double [nVar];
  }
  
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjTNE2();
  
}

CCentLax_AdjTNE2::~CCentLax_AdjTNE2(void) {
  
  delete [] Normal_ij;
  delete [] Normal_ji;
	delete [] DiffPsi;
  delete [] MeanPsi;
  delete [] MeanPsiRho;
  delete [] MeanPhi;
	delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Density_i;
  delete [] Density_j;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Proj_Jac_Tensor_i[iVar];
    delete [] Proj_Jac_Tensor_j[iVar];
  }
  delete [] Proj_Jac_Tensor_i;
  delete [] Proj_Jac_Tensor_j;
}

void CCentLax_AdjTNE2::ComputeResidual (double *val_resconv_i,
                                        double *val_resvisc_i,
                                        double *val_resconv_j,
                                        double *val_resvisc_j,
                                        double **val_Jacobian_ii,
                                        double **val_Jacobian_ij,
                                        double **val_Jacobian_ji,
                                        double **val_Jacobian_jj,
                                        CConfig *config) {
  bool ionization;
  unsigned short iDim, jDim, iSpecies, iVar, jVar, nHeavy, nEl;
  double Energy_ve_i, Energy_ve_j, rhoCvtr_i, rhoCvtr_j, rhoCvve_i, rhoCvve_j;
  double *Ms, Ru, conc_i, conc_j, rho_el_i, rho_el_j;
  double dPdrhoE_i, dPdrhoE_j, dPdrhoEve_i, dPdrhoEve_j;
  
  /*--- Initialize the residuals ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_resconv_i[iVar] = 0.0;
    val_resconv_i[iVar] = 0.0;
  }
  
  /*--- Read parameters from config ---*/
  ionization = config->GetIonization();
  
  /*--- Compute face area ---*/
	Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
  /*--- Compute and unit normal vector ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		UnitNormal[iDim] = Normal[iDim]/Area;
    Normal_ij[iDim] = Normal[iDim];
    Normal_ji[iDim] = -Normal[iDim];
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Calculate the mean values of the adjoint variables ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    MeanPsi[iVar] = 0.5 * (Psi_i[iVar] + Psi_j[iVar]);
    DiffPsi[iVar] = Psi_i[iVar]-Psi_j[iVar];
  }
  
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    MeanPsiRho[iSpecies] =  0.5*(Psi_i[iSpecies]+Psi_j[iSpecies]);
//	for (iDim = 0; iDim < nDim; iDim++)
//		MeanPhi[iDim] = 0.5 * (Psi_i[nSpecies+iDim]   + Psi_j[nSpecies+iDim]);
//	MeanPsiE        = 0.5 * (Psi_i[nSpecies+nDim]   + Psi_j[nSpecies+nDim]);
//  MeanPsiEve      = 0.5 * (Psi_i[nSpecies+nDim+1] + Psi_j[nSpecies+nDim+1]);
  
  /*--- Determine the number of heavy particle species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    Density_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = V_i[VEL_INDEX+iDim];
		Velocity_j[iDim] = V_j[VEL_INDEX+iDim];
	}
	Pressure_i   = V_i[P_INDEX];    Pressure_j   = V_j[P_INDEX];
  Enthalpy_i   = V_i[H_INDEX];    Enthalpy_j   = V_j[H_INDEX];
  SoundSpeed_i = V_i[A_INDEX];    SoundSpeed_j = V_j[A_INDEX];
  Energy_ve_i  = U_i[nSpecies+nDim+1] / V_i[RHO_INDEX];
  Energy_ve_j  = U_j[nSpecies+nDim+1] / V_j[RHO_INDEX];
  rhoCvtr_i    = V_i[RHOCVTR_INDEX];
  rhoCvtr_j    = V_j[RHOCVTR_INDEX];
  rhoCvve_i    = V_i[RHOCVVE_INDEX];
  rhoCvve_j    = V_j[RHOCVVE_INDEX];
  
  /*--- Read parameters from CConfig ---*/
  Ms = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  
  /*--- Calculate useful quantities ---*/
  if (ionization) {
    rho_el_i = Density_i[nSpecies-1];
    rho_el_j = Density_j[nSpecies-1];
  } else {
    rho_el_i = 0.0;
    rho_el_j = 0.0;
  }
  conc_i = 0.0;
  conc_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc_i += Density_i[iSpecies] / Ms[iSpecies];
    conc_j += Density_j[iSpecies] / Ms[iSpecies];
  }
  dPdrhoE_i   = conc_i*Ru / rhoCvtr_i;
  dPdrhoE_j   = conc_j*Ru / rhoCvtr_j;
  dPdrhoEve_i = -dPdrhoE_i + rho_el_i * Ru/Ms[nSpecies-1] * 1.0/rhoCvve_i;
  dPdrhoEve_j = -dPdrhoE_j + rho_el_j * Ru/Ms[nSpecies-1] * 1.0/rhoCvve_j;
  
	/*--- Calculate Projected Flux Jacobians (inviscid) ---*/
  GetInviscidProjJac(Density_i, Velocity_i, &Enthalpy_i, &Energy_ve_i,
                     dPdrhos_i, dPdrhoE_i, dPdrhoEve_i, Normal_ij, 1.0,
                     Proj_Jac_Tensor_i);
  GetInviscidProjJac(Density_j, Velocity_j, &Enthalpy_j, &Energy_ve_j,
                     dPdrhos_j, dPdrhoE_j, dPdrhoEve_j, Normal_ji, 1.0,
                     Proj_Jac_Tensor_j);
	
  
  /*--- Compute inviscid residual at point i ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_resconv_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*MeanPsi[jVar];
      val_resconv_j[iVar] += Proj_Jac_Tensor_j[jVar][iVar]*MeanPsi[jVar];
    }
  }
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_ij[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_jj[iVar][jVar] = Proj_Jac_Tensor_j[jVar][iVar];
        val_Jacobian_ji[iVar][jVar] = Proj_Jac_Tensor_j[jVar][iVar];
      }
    }
  }
  
	/*--- Compute spectral radius ---*/
  ProjVelocity_i = 0.0;
  ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
  }
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	/*--- Compute streching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc2*double(nDim)/3.0;
  
	/*--- Artifical dissipation evaluation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = Epsilon_0*StretchingFactor*MeanLambda*DiffPsi[iVar];
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
	}
  
	/*--- Contribution to implicit part ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
		}
	}
  
}
