/*!
 * \file numerics_direct_tne2.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.0.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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
#include "../include/variable_structure.hpp"
#include <limits>


CUpwRoe_TNE2::CUpwRoe_TNE2(unsigned short val_nDim, unsigned short val_nVar,
                           unsigned short val_nPrimVar,
                           unsigned short val_nPrimVarGrad,
                           CConfig *config) : CNumerics(val_nDim, val_nVar,
                                                        config) {
	
	unsigned short iVar;
  
  /*--- Read configuration parameters ---*/
	implicit   = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  ionization = config->GetIonization();
	//kappa      =config->GetRoe_Kappa()    
	
  /*--- Define useful constants ---*/
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  
  /*--- Allocate arrays ---*/
	Diff_U      = new su2double  [nVar];
  RoeU        = new su2double  [nVar];
  RoeV        = new su2double  [nPrimVar];
  RoedPdU     = new su2double  [nVar];
  RoeEve      = new su2double  [nSpecies];
	Lambda      = new su2double  [nVar];
	Epsilon     = new su2double  [nVar];
	P_Tensor    = new su2double* [nVar];
	invP_Tensor = new su2double* [nVar];
	
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]    = new su2double [nVar];
		invP_Tensor[iVar] = new su2double [nVar];
	}
  
	ProjFlux_i = new su2double [nVar];
	ProjFlux_j = new su2double [nVar];
  
}

CUpwRoe_TNE2::~CUpwRoe_TNE2(void) {
	
	unsigned short iVar;
  
	delete [] Diff_U;
  delete [] RoeU;
  delete [] RoeV;
  delete [] RoedPdU;
  delete [] RoeEve;
	delete [] Lambda;
	delete [] Epsilon;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
	}
	
	delete [] P_Tensor;
	delete [] invP_Tensor;
  delete [] ProjFlux_i;
	delete [] ProjFlux_j;
}

void CUpwRoe_TNE2::ComputeResidual(su2double *val_residual,
                                   su2double **val_Jacobian_i,
                                   su2double **val_Jacobian_j,
                                   CConfig *config) {
																	  
  unsigned short iDim, iSpecies, iVar, jVar, kVar;
  
  /*--- Compute geometrical quantities ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
	
	for (iDim = 0; iDim < nDim; iDim++)
		UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Calculate Roe variables ---*/
  R    = sqrt(abs(V_j[RHO_INDEX]/V_i[RHO_INDEX]));
  
	for (iVar = 0; iVar < nVar; iVar++)
    RoeU[iVar] = (R*U_j[iVar] + U_i[iVar])/(R+1);
  
	for (iVar = 0; iVar < nPrimVar; iVar++)
    RoeV[iVar] = (R*V_j[iVar] + V_i[iVar])/(R+1);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    RoeEve[iSpecies] = var->CalcEve(config, RoeV[TVE_INDEX], iSpecies);
  
  /*--- Calculate derivatives of pressure ---*/
  var->CalcdPdU(RoeV, RoeEve, config, RoedPdU);
  
  /*--- Calculate dual grid tangent vectors for P & invP ---*/
  CreateBasis(UnitNormal);
  
  /*--- Compute the inviscid projected fluxes ---*/
  GetInviscidProjFlux(U_i, V_i, Normal, ProjFlux_i);
  GetInviscidProjFlux(U_j, V_j, Normal, ProjFlux_j);
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix(RoeU, RoeV, RoedPdU, UnitNormal, l, m, P_Tensor);
  GetPMatrix_inv(RoeU, RoeV, RoedPdU, UnitNormal, l, m, invP_Tensor);
  
  /*--- Compute projected velocities ---*/
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeV[VEL_INDEX+iDim] * UnitNormal[iDim];
    ProjVelocity_i += V_i[VEL_INDEX+iDim]  * UnitNormal[iDim];
    ProjVelocity_j += V_j[VEL_INDEX+iDim]  * UnitNormal[iDim];
  }
  
  RoeSoundSpeed = sqrt((1.0+RoedPdU[nSpecies+nDim])*
                       RoeV[P_INDEX]/RoeV[RHO_INDEX]);
  
  /*--- Calculate eigenvalues ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Lambda[iSpecies] = ProjVelocity;
  
	for (iDim = 0; iDim < nDim-1; iDim++)
    Lambda[nSpecies+iDim] = ProjVelocity;
  
	Lambda[nSpecies+nDim-1] = ProjVelocity + RoeSoundSpeed;
  Lambda[nSpecies+nDim]   = ProjVelocity - RoeSoundSpeed;
  Lambda[nSpecies+nDim+1] = ProjVelocity;
  
  /*--- Harten and Hyman (1983) entropy correction ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Epsilon[iSpecies] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
                                         ProjVelocity_j-Lambda[iDim] ));
  for (iDim = 0; iDim < nDim-1; iDim++)
    Epsilon[nSpecies+iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
                                              ProjVelocity_j-Lambda[iDim] ));
  Epsilon[nSpecies+nDim-1] = 4.0*max(0.0, max(Lambda[nSpecies+nDim-1]-(ProjVelocity_i+V_i[A_INDEX]),
                                              (ProjVelocity_j+V_j[A_INDEX])-Lambda[nSpecies+nDim-1]));
  Epsilon[nSpecies+nDim]   = 4.0*max(0.0, max(Lambda[nSpecies+nDim]-(ProjVelocity_i-V_i[A_INDEX]),
                                              (ProjVelocity_j-V_j[A_INDEX])-Lambda[nSpecies+nDim]));
  Epsilon[nSpecies+nDim+1] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
                                              ProjVelocity_j-Lambda[iDim] ));
  for (iVar = 0; iVar < nVar; iVar++)
    if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
      Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
    else
      Lambda[iVar] = fabs(Lambda[iVar]);
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  /*--- Calculate inviscid projected Jacobians ---*/
  // Note: Scaling value is 0.5 because inviscid flux is based on 0.5*(Fc_i+Fc_j)
  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 0.5, val_Jacobian_i);
  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 0.5, val_Jacobian_j);
    
  /*--- Difference of conserved variables at iPoint and jPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
  
  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5 * (ProjFlux_i[iVar] + ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {

      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
  
	    val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
      val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
    }
  }
}