/*!
 * \file ausm.cpp
 * \brief Implementations of the AUSM-family of schemes in NEMO.
 * \author F. Palacios, S.R. Copeland, W. Maier, C. Garbacz
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/NEMO/convection/ausm.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSM_NEMO::CUpwAUSM_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad,
                             CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  FcL    = new su2double [nVar];
  FcR    = new su2double [nVar];
  dmLP   = new su2double [nVar];
  dmRM   = new su2double [nVar];
  dpLP   = new su2double [nVar];
  dpRM   = new su2double [nVar];
  rhos_i = new su2double [nSpecies];
  rhos_j = new su2double [nSpecies];
  u_i    = new su2double [nDim];
  u_j    = new su2double [nDim];
  daL    = new su2double [nVar];
  daR    = new su2double [nVar];

  Flux   = new su2double[nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwAUSM_NEMO::~CUpwAUSM_NEMO(void) {

  delete [] FcL;
  delete [] FcR;
  delete [] dmLP;
  delete [] dmRM;
  delete [] dpLP;
  delete [] dpRM;
  delete [] rhos_i;
  delete [] rhos_j;
  delete [] u_i;
  delete [] u_j;
  delete [] Flux;
  delete [] daL;
  delete [] daR;
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
}

CNumerics::ResidualType<> CUpwAUSM_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iDim, iVar, jVar, iSpecies;
  su2double rho_i, rho_j, 
  e_ve_i, e_ve_j, mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;

  /*--- Compute geometric quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim] = V_i[VEL_INDEX+iDim];
    u_j[iDim] = V_j[VEL_INDEX+iDim];
  }

  P_i   = V_i[P_INDEX];   P_j   = V_j[P_INDEX];
  h_i   = V_i[H_INDEX];   h_j   = V_j[H_INDEX];
  a_i   = V_i[A_INDEX];   a_j   = V_j[A_INDEX];
  rho_i = V_i[RHO_INDEX]; rho_j = V_j[RHO_INDEX];

  e_ve_i  = 0; e_ve_j  = 0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    e_ve_i += (V_i[RHOS_INDEX+iSpecies]*eve_i[iSpecies])/rho_i;
    e_ve_j += (V_j[RHOS_INDEX+iSpecies]*eve_j[iSpecies])/rho_j;
  }

  /*--- Projected velocities ---*/
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }

  /*--- Calculate L/R Mach numbers ---*/
  mL = ProjVel_i/a_i;
  mR = ProjVel_j/a_j;

  /*--- Calculate split numerical fluxes ---*/
  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
  else                 mLP = 0.5*(mL+fabs(mL));

  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
  else                 mRM = 0.5*(mR-fabs(mR));

  mF = mLP + mRM;

  if (fabs(mL) <= 1.0) pLP = 0.25*P_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
  else                 pLP = 0.5*P_i*(mL+fabs(mL))/mL;

  if (fabs(mR) <= 1.0) pRM = 0.25*P_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
  else                 pRM = 0.5*P_j*(mR-fabs(mR))/mR;

  pF = pLP + pRM;
  Phi = fabs(mF);

  /*--- Assign left & right convective vectors ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    FcL[iSpecies] = rhos_i[iSpecies]*a_i;
    FcR[iSpecies] = rhos_j[iSpecies]*a_j;
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    FcL[nSpecies+iDim] = rho_i*a_i*u_i[iDim];
    FcR[nSpecies+iDim] = rho_j*a_j*u_j[iDim];
  }
  FcL[nSpecies+nDim]   = rho_i*a_i*h_i;
  FcR[nSpecies+nDim]   = rho_j*a_j*h_j;
  FcL[nSpecies+nDim+1] = rho_i*a_i*e_ve_i;
  FcR[nSpecies+nDim+1] = rho_j*a_j*e_ve_j;

  /*--- Compute numerical flux ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = 0.5*((mF+Phi)*FcL[iVar]+(mF-Phi)*FcR[iVar])*Area;

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[nSpecies+iDim] += pF*UnitNormal[iDim]*Area;

  /********************************************************/
  //TODO DELETE THIS: ROE Jacobian for DEBUGGING
  if (implicit){

    RoeU        = new su2double  [nVar];
    RoeV        = new su2double  [nPrimVar];
    RoedPdU     = new su2double  [nVar];
    Lambda      = new su2double  [nVar];
    Epsilon     = new su2double  [nVar];
    P_Tensor    = new su2double* [nVar];
    invP_Tensor = new su2double* [nVar];

    roe_eves.resize(nSpecies,0.0);

    for (iVar = 0; iVar < nVar; iVar++) {
      P_Tensor[iVar]    = new su2double [nVar];
      invP_Tensor[iVar] = new su2double [nVar];
    }
	  
    su2double R = sqrt(fabs(rho_i/rho_j));
    
    
    for (iVar = 0; iVar < nVar; iVar++){
      RoeU[iVar] = (R*U_j[iVar] + U_i[iVar])/(R+1);
    }
    for (iVar = 0; iVar < nPrimVar; iVar++){
      RoeV[iVar] = (R*V_j[iVar] + V_i[iVar])/(R+1);
    }
    
    auto& roe_eves = fluidmodel->ComputeSpeciesEve(RoeV[TVE_INDEX]);

    /*--- Calculate derivatives of pressure ---*/
    fluidmodel->ComputedPdU(RoeV, roe_eves, RoedPdU);

    /*--- Calculate dual grid tangent vectors for P & invP ---*/
    su2double l[MAXNDIM], m[MAXNDIM];
    CreateBasis(UnitNormal,l,m);

    /*--- Compute projected P, invP, and Lambda ---*/
    GetPMatrix    (RoeU, RoeV, RoedPdU, UnitNormal, l, m, P_Tensor);
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

    /*--- Calculate inviscid projected Jacobians ---*/
    // Note: Scaling value is 0.5 because inviscid flux is based on 0.5*(Fc_i+Fc_j)
    if (implicit){
      GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 0.5, Jacobian_i);
      GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 0.5, Jacobian_j);
    }

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {

        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        Proj_ModJac_Tensor_ij = 0.0;
        for (unsigned short kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

        Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        delete [] P_Tensor[iVar];
        delete [] invP_Tensor[iVar];
      }
      
      delete [] RoeU;
      delete [] RoeV;
      delete [] RoedPdU;
      delete [] Lambda;
      delete [] Epsilon;
      delete [] P_Tensor;
      delete [] invP_Tensor;
  }
  /*********************************************************/


  // if (implicit){

  //   auto& Ms   = fluidmodel->GetSpeciesMolarMass();
  //   auto& Cvtr = fluidmodel->GetSpeciesCvTraRot();
  //   Ru         = 1000.0*UNIVERSAL_GAS_CONSTANT;
  //   rhoCvtr_i  = V_i[RHOCVTR_INDEX];
  //   rhoCvtr_j  = V_j[RHOCVTR_INDEX];

  //   /*--- Initialize the Jacobians ---*/
  //   for (iVar = 0; iVar < nVar; iVar++) {
  //     for (jVar = 0; jVar < nVar; jVar++) {
  //       Jacobian_i[iVar][jVar] = 0.0;
  //       Jacobian_j[iVar][jVar] = 0.0;
  //     }
  //   }

  //   if (mF >= 0.0) FcLR = FcL;
  //   else           FcLR = FcR;

  //   /*--- Sound speed derivatives: Species density ---*/
  //   for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
  //     daL[iSpecies] = 1.0/(2.0*a_i) * (1/rhoCvtr_i*(Ru/Ms[iSpecies] - Cvtr[iSpecies]*dPdU_i[nSpecies+nDim])*P_i/rho_i
  //                   + 1.0/rho_i*(1.0+dPdU_i[nSpecies+nDim])*(dPdU_i[iSpecies] - P_i/rho_i));
  //     daR[iSpecies] = 1.0/(2.0*a_j) * (1/rhoCvtr_j*(Ru/Ms[iSpecies] - Cvtr[iSpecies]*dPdU_j[nSpecies+nDim])*P_j/rho_j
  //                   + 1.0/rho_j*(1.0+dPdU_j[nSpecies+nDim])*(dPdU_j[iSpecies] - P_j/rho_j));
  //   }
  //   for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
  //     daL[nSpecies-1] = 1.0/(2.0*a_i*rho_i) * (1+dPdU_i[nSpecies+nDim])*(dPdU_i[nSpecies-1] - P_i/rho_i);
  //     daR[nSpecies-1] = 1.0/(2.0*a_j*rho_j) * (1+dPdU_j[nSpecies+nDim])*(dPdU_j[nSpecies-1] - P_j/rho_j);
  //   }

  //   /*--- Sound speed derivatives: Momentum ---*/
  //   for (iDim = 0; iDim < nDim; iDim++) {
  //     daL[nSpecies+iDim] = -1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim])*u_i[iDim];
  //     daR[nSpecies+iDim] = -1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim])*u_j[iDim];
  //   }

  //   /*--- Sound speed derivatives: Energy ---*/
  //   daL[nSpecies+nDim]   = 1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim]);
  //   daR[nSpecies+nDim]   = 1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim]);

  //   /*--- Sound speed derivatives: Vib-el energy ---*/
  //   daL[nSpecies+nDim+1] = 1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim+1]);
  //   daR[nSpecies+nDim+1] = 1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim+1]);

  //   /*--- Left state Jacobian ---*/
  //   if (mF >= 0) {

  //     /*--- Jacobian contribution: dFc terms ---*/
  //     for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
  //       for (jVar = 0; jVar < nVar; jVar++) {
  //         Jacobian_i[iVar][jVar] += mF * FcL[iVar]/a_i * daL[jVar];
  //       }
  //       Jacobian_i[iVar][iVar] += mF * a_i;
  //     }
  //     for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //       Jacobian_i[nSpecies+nDim][iSpecies] += mF * (dPdU_i[iSpecies]*a_i + rho_i*h_i*daL[iSpecies]);
  //     }
  //     for (iDim = 0; iDim < nDim; iDim++) {
  //       Jacobian_i[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_i[nSpecies+nDim]*u_i[iDim]*a_i + rho_i*h_i*daL[nSpecies+iDim]);
  //     }
  //     Jacobian_i[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_i[nSpecies+nDim])*a_i + rho_i*h_i*daL[nSpecies+nDim]);
  //     Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_i[nSpecies+nDim+1]*a_i + rho_i*h_i*daL[nSpecies+nDim+1]);
  //     for (jVar = 0; jVar < nVar; jVar++) {
  //       Jacobian_i[nSpecies+nDim+1][jVar] +=  mF * FcL[nSpecies+nDim+1]/a_i * daL[jVar];
  //     }
  //     Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += mF * a_i;
  //   }

  //   /*--- Calculate derivatives of the split pressure flux ---*/
  //   if ( (mF >= 0) || ((mF < 0)&&(fabs(mF) <= 1.0)) ) {
  //     if (fabs(mL) <= 1.0) {

  //       /*--- Mach number ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dmLP[iSpecies] = 0.5*(mL+1.0) * (-ProjVel_i/(rho_i*a_i) - ProjVel_i*daL[iSpecies]/(a_i*a_i));
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dmLP[nSpecies+iDim] = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*a_i));
  //       dmLP[nSpecies+nDim]   = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim]);
  //       dmLP[nSpecies+nDim+1] = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1]);

  //       /*--- Pressure ---*/
  //       for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dpLP[iSpecies] = 0.25*(mL+1.0) * (dPdU_i[iSpecies]*(mL+1.0)*(2.0-mL)
  //                                         + P_i*(-ProjVel_i/(rho_i*a_i)
  //                                                -ProjVel_i*daL[iSpecies]/(a_i*a_i))*(3.0-3.0*mL));
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dpLP[nSpecies+iDim] = 0.25*(mL+1.0) * (-u_i[iDim]*dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
  //                                               + P_i*( -ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim]
  //                                               + UnitNormal[iDim]/(rho_i*a_i))*(3.0-3.0*mL));
  //       dpLP[nSpecies+nDim]   = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
  //                                               + P_i*(-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim])*(3.0-3.0*mL));
  //       dpLP[nSpecies+nDim+1] = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim+1]*(mL+1.0)*(2.0-mL)
  //                                               + P_i*(-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1])*(3.0-3.0*mL));
  //     } else {

  //       /*--- Mach number ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dmLP[iSpecies]      = -ProjVel_i/(rho_i*a_i) - ProjVel_i*daL[iSpecies]/(a_i*a_i);
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dmLP[nSpecies+iDim] = -ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*a_i);
  //       dmLP[nSpecies+nDim]   = -ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim];
  //       dmLP[nSpecies+nDim+1] = -ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1];

  //       /*--- Pressure ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dpLP[iSpecies] = dPdU_i[iSpecies];
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dpLP[nSpecies+iDim] = (-u_i[iDim]*dPdU_i[nSpecies+nDim]);
  //       dpLP[nSpecies+nDim]   = dPdU_i[nSpecies+nDim];
  //       dpLP[nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1];
  //     }

  //     /*--- dM contribution ---*/
  //     for (iVar = 0; iVar < nVar; iVar++) {
  //       for (jVar = 0; jVar < nVar; jVar++) {
  //         Jacobian_i[iVar][jVar] += dmLP[jVar]*FcLR[iVar];
  //       }
  //     }

  //     /*--- Jacobian contribution: dP terms ---*/
  //     for (iDim = 0; iDim < nDim; iDim++) {
  //       for (iVar = 0; iVar < nVar; iVar++) {
  //         Jacobian_i[nSpecies+iDim][iVar] += dpLP[iVar]*UnitNormal[iDim];
  //       }
  //     }
  //   }

  //   /*--- Right state Jacobian ---*/
  //   if (mF < 0) {

  //     /*--- Jacobian contribution: dFc terms ---*/
  //     for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
  //       for (jVar = 0; jVar < nVar; jVar++) {
  //         Jacobian_j[iVar][jVar] += mF * FcR[iVar]/a_j * daR[jVar];
  //       }
  //       Jacobian_j[iVar][iVar] += mF * a_j;
  //     }
  //     for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //       Jacobian_j[nSpecies+nDim][iSpecies] += mF * (dPdU_j[iSpecies]*a_j + rho_j*h_j*daR[iSpecies]);
  //     }
  //     for (iDim = 0; iDim < nDim; iDim++) {
  //       Jacobian_j[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_j[nSpecies+nDim]*u_j[iDim]*a_j + rho_j*h_j*daR[nSpecies+iDim]);
  //     }
  //     Jacobian_j[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_j[nSpecies+nDim])*a_j + rho_j*h_j*daR[nSpecies+nDim]);
  //     Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_j[nSpecies+nDim+1]*a_j + rho_j*h_j*daR[nSpecies+nDim+1]);
  //     for (jVar = 0; jVar < nVar; jVar++) {
  //       Jacobian_j[nSpecies+nDim+1][jVar] +=  mF * FcR[nSpecies+nDim+1]/a_j * daR[jVar];
  //     }
  //     Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += mF * a_j;
  //   }

  //   /*--- Calculate derivatives of the split pressure flux ---*/
  //   if ( (mF < 0) || ((mF >= 0)&&(fabs(mF) <= 1.0)) ) {
  //     if (fabs(mR) <= 1.0) {

  //       /*--- Mach ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dmRM[iSpecies] = -0.5*(mR-1.0) * (-ProjVel_j/(rho_j*a_j) - ProjVel_j*daR[iSpecies]/(a_j*a_j));
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dmRM[nSpecies+iDim] = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*a_j));
  //       dmRM[nSpecies+nDim]   = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim]);
  //       dmRM[nSpecies+nDim+1] = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1]);

  //       /*--- Pressure ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dpRM[iSpecies] = 0.25*(mR-1.0) * (dPdU_j[iSpecies]*(mR-1.0)*(2.0+mR)
  //                                         + P_j*(-ProjVel_j/(rho_j*a_j)
  //                                                -ProjVel_j*daR[iSpecies]/(a_j*a_j))*(3.0+3.0*mR));
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dpRM[nSpecies+iDim] = 0.25*(mR-1.0) * ((-u_j[iDim]*dPdU_j[nSpecies+nDim])*(mR-1.0)*(2.0+mR)
  //                                                + P_j*( -ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim]
  //                                                + UnitNormal[iDim]/(rho_j*a_j))*(3.0+3.0*mR));
  //       dpRM[nSpecies+nDim]   = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim]*(mR-1.0)*(2.0+mR)
  //                                              + P_j*(-ProjVel_j/(a_j*a_j)*daR[nSpecies+nDim])*(3.0+3.0*mR));
  //       dpRM[nSpecies+nDim+1] = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim+1]*(mR-1.0)*(2.0+mR)
  //                                              + P_j*(-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1])*(3.0+3.0*mR));

  //     } else {

  //       /*--- Mach ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dmRM[iSpecies]      = -ProjVel_j/(rho_j*a_j) - ProjVel_j*daR[iSpecies]/(a_j*a_j);
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dmRM[nSpecies+iDim] = -ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*a_j);
  //       dmRM[nSpecies+nDim]   = -ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim];
  //       dmRM[nSpecies+nDim+1] = -ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1];

  //       /*--- Pressure ---*/
  //       for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //         dpRM[iSpecies] = dPdU_j[iSpecies];
  //       for (iDim = 0; iDim < nDim; iDim++)
  //         dpRM[nSpecies+iDim] = -u_j[iDim]*dPdU_j[nSpecies+nDim];
  //       dpRM[nSpecies+nDim]   = dPdU_j[nSpecies+nDim];
  //       dpRM[nSpecies+nDim+1] = dPdU_j[nSpecies+nDim+1];
  //     }

  //     /*--- Jacobian contribution: dM terms ---*/
  //     for (iVar = 0; iVar < nVar; iVar++) {
  //       for (jVar = 0; jVar < nVar; jVar++) {
  //         Jacobian_j[iVar][jVar] += dmRM[jVar] * FcLR[iVar];
  //       }
  //     }

  //     /*--- Jacobian contribution: dP terms ---*/
  //     for (iDim = 0; iDim < nDim; iDim++) {
  //       for (iVar = 0; iVar < nVar; iVar++) {
  //         Jacobian_j[nSpecies+iDim][iVar] += dpRM[iVar]*UnitNormal[iDim];
  //       }
  //     }
  //   }

  //   /*--- Integrate over dual-face area ---*/
  //   for (iVar = 0; iVar < nVar; iVar++) {
  //     for (jVar = 0; jVar < nVar; jVar++) {
  //       Jacobian_i[iVar][jVar] *= Area;
  //       Jacobian_j[iVar][jVar] *= Area;
  //     }
  //   }
  // }
  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}
