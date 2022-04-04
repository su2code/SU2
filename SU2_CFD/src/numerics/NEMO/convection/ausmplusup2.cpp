/*!
 * \file ausmplusup2.cpp
 * \brief Implementations of the AUSM-family of schemes - AUSM+UP2.
 * \author W. Maier, A. Sachedeva, C. Garbacz
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/NEMO/convection/ausmplusup2.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSMPLUSUP2_NEMO::CUpwAUSMPLUSUP2_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                                           unsigned short val_nPrimVar,
                                           unsigned short val_nPrimVarGrad,
                                           CConfig *config): CNEMONumerics (val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config){

  /*--- Define useful constants ---*/
  Kp       = 0.25;
  sigma    = 1.0;

  /*--- Allocate data structures ---*/
  FcL    = new su2double [nVar];
  FcR    = new su2double [nVar];
  dmLP   = new su2double [nVar];
  dmRM   = new su2double [nVar];
  dpLP   = new su2double [nVar];
  dpRM   = new su2double [nVar];
  daL    = new su2double [nVar];
  daR    = new su2double [nVar];
  rhos_i = new su2double [nSpecies];
  rhos_j = new su2double [nSpecies];
  u_i    = new su2double [nDim];
  u_j    = new su2double [nDim];

  Flux   = new su2double[nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwAUSMPLUSUP2_NEMO::~CUpwAUSMPLUSUP2_NEMO(void) {

  delete [] FcL;
  delete [] FcR;
  delete [] dmLP;
  delete [] dmRM;
  delete [] dpLP;
  delete [] dpRM;
  delete [] daL;
  delete [] daR;
  delete [] rhos_i;
  delete [] rhos_j;
  delete [] u_i;
  delete [] u_j;
  delete [] Flux;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
}

CNumerics::ResidualType<> CUpwAUSMPLUSUP2_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iDim, iVar,jVar, iSpecies;
  su2double rho_i, rho_j,
  mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi,
  e_ve_i, e_ve_j;

  /*--- Face area ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  Minf  = config->GetMach();

  /*--- Extracting primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim]  = V_i[VEL_INDEX+iDim];
    u_j[iDim]  = V_j[VEL_INDEX+iDim];
  }
  su2double sq_veli = GeometryToolbox::SquaredNorm(nDim, u_i);
  su2double sq_velj = GeometryToolbox::SquaredNorm(nDim, u_j);

  P_i   = V_i[P_INDEX];   P_j   = V_j[P_INDEX];
  h_i   = V_i[H_INDEX];   h_j   = V_j[H_INDEX];
  a_i   = V_i[A_INDEX];   a_j   = V_j[A_INDEX];
  rho_i = V_i[RHO_INDEX]; rho_j = V_j[RHO_INDEX];

  rhoCvtr_i = V_i[RHOCVTR_INDEX]; rhoCvtr_j = V_j[RHOCVTR_INDEX];
  rhoCvve_i = V_i[RHOCVVE_INDEX]; rhoCvve_j = V_j[RHOCVVE_INDEX];

  e_ve_i = 0.0; e_ve_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    e_ve_i += (rhos_i[iSpecies]*eve_i[iSpecies])/rho_i;
    e_ve_j += (rhos_j[iSpecies]*eve_j[iSpecies])/rho_j;
  }

  /*--- Projected velocities ---*/
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }

  /*--- Compute C*  ---*/
  CstarL = sqrt(2.0*(Gamma_i-1.0)/(Gamma_i+1.0)*h_i);
  CstarR = sqrt(2.0*(Gamma_j-1.0)/(Gamma_j+1.0)*h_j);

  /*--- Compute C^ ---*/
  ChatL = CstarL*CstarL/max(CstarL,ProjVel_i);
  ChatR = CstarR*CstarR/max(CstarR,-ProjVel_j);

  /*--- Interface speed of sound ---*/
  aF = min(ChatL,ChatR);

  mL  = ProjVel_i/aF;
  mR  = ProjVel_j/aF;

  rhoF = 0.5*(rho_i+rho_j);
  MFsq = 0.5*(mL*mL+mR*mR);

  param1 = max(MFsq, Minf*Minf);
  Mrefsq = (min(1.0, param1));
  fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  beta = 1.0/8.0;

  /*--- Pressure diffusion term ---*/
  Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(P_j-P_i)/(rhoF*aF*aF);

  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0)+beta*(mL*mL-1.0)*(mL*mL-1.0);
  else                 mLP = 0.5*(mL+fabs(mL));

  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0)-beta*(mR*mR-1.0)*(mR*mR-1.0);
  else                 mRM = 0.5*(mR-fabs(mR));

  mF = mLP + mRM + Mp;

  if (fabs(mL) <= 1.0) pLP = (0.25*(mL+1.0)*(mL+1.0)*(2.0-mL)+alpha*mL*(mL*mL-1.0)*(mL*mL-1.0));
  else                 pLP = 0.5*(mL+fabs(mL))/mL;

  if (fabs(mR) <= 1.0) pRM = (0.25*(mR-1.0)*(mR-1.0)*(2.0+mR)-alpha*mR*(mR*mR-1.0)*(mR*mR-1.0));
  else                 pRM = 0.5*(mR-fabs(mR))/mR;

  /*... Modified pressure flux ...*/
  //Use this definition
  pFi = sqrt(0.5*(sq_veli+sq_velj))*(pLP+pRM-1.0)*0.5*(rho_j+rho_i)*aF;
  pF  = 0.5*(P_j+P_i)+0.5*(pLP-pRM)*(P_i-P_j)+pFi;

  Phi = fabs(mF);

  mfP=0.5*(mF+Phi);
  mfM=0.5*(mF-Phi);

  /*--- Assign left & right covective fluxes ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    FcL[iSpecies] = rhos_i[iSpecies]*aF;
    FcR[iSpecies] = rhos_j[iSpecies]*aF;
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    FcL[nSpecies+iDim] = rho_i*aF*u_i[iDim];
    FcR[nSpecies+iDim] = rho_j*aF*u_j[iDim];
  }
  FcL[nSpecies+nDim]   = rho_i*aF*h_i;
  FcR[nSpecies+nDim]   = rho_j*aF*h_j;
  FcL[nSpecies+nDim+1] = rho_i*aF*e_ve_i;
  FcR[nSpecies+nDim+1] = rho_j*aF*e_ve_j;

  /*--- Compute numerical flux ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = (mfP*FcL[iVar]+mfM*FcR[iVar])*Area;

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[nSpecies+iDim] += pF*UnitNormal[iDim]*Area;

  /*--- AUSM's Jacobian ---*/
  if (implicit)  {

    /*--- Initialize the Jacobians ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }
    }

    if (mF >= 0.0) FcLR = FcL;
    else           FcLR = FcR;

    /*--- Sound speed derivatives: Species density ---*/
    const auto& Cvtrs = fluidmodel -> GetSpeciesCvTraRot();
    const auto& Ms    = fluidmodel -> GetSpeciesMolarMass();
    su2double Ru = 1000.0*UNIVERSAL_GAS_CONSTANT;

    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
      daL[iSpecies] = 1.0/(2.0*aF*rho_i) * (1+dPdU_i[nSpecies+nDim])*(dPdU_i[iSpecies] - P_i/rho_i);
      daR[iSpecies] = 1.0/(2.0*aF*rho_j) * (1+dPdU_j[nSpecies+nDim])*(dPdU_j[iSpecies] - P_j/rho_j);
    }

    for (iSpecies = nEl; iSpecies < nHeavy; iSpecies++) {
      daL[iSpecies] = 1.0/(2.0*aF) * (1/rhoCvtr_i*(Ru/Ms[iSpecies] - Cvtrs[iSpecies]*dPdU_i[nSpecies+nDim])*P_i/rho_i
          + 1.0/rho_i*(1.0+dPdU_i[nSpecies+nDim])*(dPdU_i[iSpecies] - P_i/rho_i));
      daR[iSpecies] = 1.0/(2.0*aF) * (1/rhoCvtr_j*(Ru/Ms[iSpecies] - Cvtrs[iSpecies]*dPdU_j[nSpecies+nDim])*P_j/rho_j
          + 1.0/rho_j*(1.0+dPdU_j[nSpecies+nDim])*(dPdU_j[iSpecies] - P_j/rho_j));
    }

    /*--- Sound speed derivatives: Momentum ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      daL[nSpecies+iDim] = -1.0/(2.0*rho_i*aF) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim])*u_i[iDim];
      daR[nSpecies+iDim] = -1.0/(2.0*rho_j*aF) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim])*u_j[iDim];
    }

    /*--- Sound speed derivatives: Energy ---*/
    daL[nSpecies+nDim]   = 1.0/(2.0*rho_i*aF) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim]);
    daR[nSpecies+nDim]   = 1.0/(2.0*rho_j*aF) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim]);

    /*--- Sound speed derivatives: Vib-el energy ---*/
    daL[nSpecies+nDim+1] = 1.0/(2.0*rho_i*aF) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim+1]);
    daR[nSpecies+nDim+1] = 1.0/(2.0*rho_j*aF) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim+1]);

    /*--- Left state Jacobian ---*/
    if (mF >= 0) {

      /*--- Jacobian contribution: dFc terms ---*/
      for (auto iVar = 0u; iVar < nSpecies+nDim; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] += mF * FcL[iVar]/aF * daL[jVar];
        }
        Jacobian_i[iVar][iVar] += mF * aF;
      }
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Jacobian_i[nSpecies+nDim][iSpecies] += mF * (dPdU_i[iSpecies]*aF + rho_i*h_i*daL[iSpecies]);
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Jacobian_i[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_i[nSpecies+nDim]*u_i[iDim]*aF + rho_i*h_i*daL[nSpecies+iDim]);
      }
      Jacobian_i[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_i[nSpecies+nDim])*aF + rho_i*h_i*daL[nSpecies+nDim]);
      Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_i[nSpecies+nDim+1]*aF + rho_i*h_i*daL[nSpecies+nDim+1]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[nSpecies+nDim+1][jVar] +=  mF * FcL[nSpecies+nDim+1]/aF * daL[jVar];
      }
      Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += mF * aF;
    }

    /*--- Calculate derivatives of the split pressure flux ---*/
    if ( (mF >= 0) || ((mF < 0)&&(fabs(mF) <= 1.0)) ) {
      if (fabs(mL) <= 1.0) {

        /*--- Mach number ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dmLP[iSpecies] = 0.5*(mL+1.0) * (-ProjVel_i/(rho_i*aF) - ProjVel_i*daL[iSpecies]/(aF*aF));
        for (iDim = 0; iDim < nDim; iDim++)
          dmLP[nSpecies+iDim] = 0.5*(mL+1.0) * (-ProjVel_i/(aF*aF) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*aF));
        dmLP[nSpecies+nDim]   = 0.5*(mL+1.0) * (-ProjVel_i/(aF*aF) * daL[nSpecies+nDim]);
        dmLP[nSpecies+nDim+1] = 0.5*(mL+1.0) * (-ProjVel_i/(aF*aF) * daL[nSpecies+nDim+1]);

        /*--- Pressure ---*/
        for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dpLP[iSpecies] = 0.25*(mL+1.0) * (dPdU_i[iSpecies]*(mL+1.0)*(2.0-mL)
                                            + P_i*(-ProjVel_i/(rho_i*aF)
                                                   -ProjVel_i*daL[iSpecies]/(aF*aF))*(3.0-3.0*mL));
        for (iDim = 0; iDim < nDim; iDim++)
          dpLP[nSpecies+iDim] = 0.25*(mL+1.0) * (-u_i[iDim]*dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
              + P_i*( -ProjVel_i/(aF*aF) * daL[nSpecies+iDim]
              + UnitNormal[iDim]/(rho_i*aF))*(3.0-3.0*mL));
        dpLP[nSpecies+nDim]   = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
            + P_i*(-ProjVel_i/(aF*aF) * daL[nSpecies+nDim])*(3.0-3.0*mL));
        dpLP[nSpecies+nDim+1] = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim+1]*(mL+1.0)*(2.0-mL)
            + P_i*(-ProjVel_i/(aF*aF) * daL[nSpecies+nDim+1])*(3.0-3.0*mL));
      } else {

        /*--- Mach number ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dmLP[iSpecies]      = -ProjVel_i/(rho_i*aF) - ProjVel_i*daL[iSpecies]/(aF*aF);
        for (iDim = 0; iDim < nDim; iDim++)
          dmLP[nSpecies+iDim] = -ProjVel_i/(aF*aF) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*aF);
        dmLP[nSpecies+nDim]   = -ProjVel_i/(aF*aF) * daL[nSpecies+nDim];
        dmLP[nSpecies+nDim+1] = -ProjVel_i/(aF*aF) * daL[nSpecies+nDim+1];

        /*--- Pressure ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dpLP[iSpecies] = dPdU_i[iSpecies];
        for (iDim = 0; iDim < nDim; iDim++)
          dpLP[nSpecies+iDim] = (-u_i[iDim]*dPdU_i[nSpecies+nDim]);
        dpLP[nSpecies+nDim]   = dPdU_i[nSpecies+nDim];
        dpLP[nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1];
      }

      /*--- dM contribution ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] += dmLP[jVar]*FcLR[iVar];
        }
      }

      /*--- Jacobian contribution: dP terms ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_i[nSpecies+iDim][iVar] += dpLP[iVar]*UnitNormal[iDim];
        }
      }
    }

    /*--- Right state Jacobian ---*/
    if (mF < 0) {

      /*--- Jacobian contribution: dFc terms ---*/
      for (auto iVar = 0u; iVar < nSpecies+nDim; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_j[iVar][jVar] += mF * FcR[iVar]/aF * daR[jVar];
        }
        Jacobian_j[iVar][iVar] += mF * aF;
      }
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Jacobian_j[nSpecies+nDim][iSpecies] += mF * (dPdU_j[iSpecies]*aF + rho_j*h_j*daR[iSpecies]);
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Jacobian_j[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_j[nSpecies+nDim]*u_j[iDim]*aF + rho_j*h_j*daR[nSpecies+iDim]);
      }
      Jacobian_j[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_j[nSpecies+nDim])*aF + rho_j*h_j*daR[nSpecies+nDim]);
      Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_j[nSpecies+nDim+1]*aF + rho_j*h_j*daR[nSpecies+nDim+1]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_j[nSpecies+nDim+1][jVar] +=  mF * FcR[nSpecies+nDim+1]/aF * daR[jVar];
      }
      Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += mF * aF;
    }

    /*--- Calculate derivatives of the split pressure flux ---*/
    if ( (mF < 0) || ((mF >= 0)&&(fabs(mF) <= 1.0)) ) {
      if (fabs(mR) <= 1.0) {

        /*--- Mach ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dmRM[iSpecies] = -0.5*(mR-1.0) * (-ProjVel_j/(rho_j*aF) - ProjVel_j*daR[iSpecies]/(aF*aF));
        for (iDim = 0; iDim < nDim; iDim++)
          dmRM[nSpecies+iDim] = -0.5*(mR-1.0) * (-ProjVel_j/(aF*aF) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*aF));
        dmRM[nSpecies+nDim]   = -0.5*(mR-1.0) * (-ProjVel_j/(aF*aF) * daR[nSpecies+nDim]);
        dmRM[nSpecies+nDim+1] = -0.5*(mR-1.0) * (-ProjVel_j/(aF*aF) * daR[nSpecies+nDim+1]);

        /*--- Pressure ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dpRM[iSpecies] = 0.25*(mR-1.0) * (dPdU_j[iSpecies]*(mR-1.0)*(2.0+mR)
                                            + P_j*(-ProjVel_j/(rho_j*aF)
                                                   -ProjVel_j*daR[iSpecies]/(aF*aF))*(3.0+3.0*mR));
        for (iDim = 0; iDim < nDim; iDim++)
          dpRM[nSpecies+iDim] = 0.25*(mR-1.0) * ((-u_j[iDim]*dPdU_j[nSpecies+nDim])*(mR-1.0)*(2.0+mR)
              + P_j*( -ProjVel_j/(aF*aF) * daR[nSpecies+iDim]
              + UnitNormal[iDim]/(rho_j*aF))*(3.0+3.0*mR));
        dpRM[nSpecies+nDim]   = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim]*(mR-1.0)*(2.0+mR)
            + P_j*(-ProjVel_j/(aF*aF)*daR[nSpecies+nDim])*(3.0+3.0*mR));
        dpRM[nSpecies+nDim+1] = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim+1]*(mR-1.0)*(2.0+mR)
            + P_j*(-ProjVel_j/(aF*aF) * daR[nSpecies+nDim+1])*(3.0+3.0*mR));

      } else {

        /*--- Mach ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dmRM[iSpecies]      = -ProjVel_j/(rho_j*aF) - ProjVel_j*daR[iSpecies]/(aF*aF);
        for (iDim = 0; iDim < nDim; iDim++)
          dmRM[nSpecies+iDim] = -ProjVel_j/(aF*aF) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*aF);
        dmRM[nSpecies+nDim]   = -ProjVel_j/(aF*aF) * daR[nSpecies+nDim];
        dmRM[nSpecies+nDim+1] = -ProjVel_j/(aF*aF) * daR[nSpecies+nDim+1];

        /*--- Pressure ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          dpRM[iSpecies] = dPdU_j[iSpecies];
        for (iDim = 0; iDim < nDim; iDim++)
          dpRM[nSpecies+iDim] = -u_j[iDim]*dPdU_j[nSpecies+nDim];
        dpRM[nSpecies+nDim]   = dPdU_j[nSpecies+nDim];
        dpRM[nSpecies+nDim+1] = dPdU_j[nSpecies+nDim+1];
      }

      /*--- Jacobian contribution: dM terms ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_j[iVar][jVar] += dmRM[jVar] * FcLR[iVar];
        }
      }

      /*--- Jacobian contribution: dP terms ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_j[nSpecies+iDim][iVar] += dpRM[iVar]*UnitNormal[iDim];
        }
      }
    }

    /*--- Integrate over dual-face area ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] *= Area;
        Jacobian_j[iVar][jVar] *= Area;
      }
    }
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}
