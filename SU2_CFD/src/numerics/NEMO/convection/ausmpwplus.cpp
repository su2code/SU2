/*!
 * \file ausmpwplus.cpp
 * \brief Implementations of the AUSM-family of schemes - AUSMPWPLUS.
 * \author F. Palacios, W. Maier, C. Garbacz
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

#include "../../../../include/numerics/NEMO/convection/ausmpwplus.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSMPWplus_NEMO::CUpwAUSMPWplus_NEMO(unsigned short val_nDim,
                                         unsigned short val_nVar,
                                         unsigned short val_nPrimVar,
                                         unsigned short val_nPrimVarGrad,
                                         CConfig *config) : CNEMONumerics(val_nDim, val_nVar,
                                                                          val_nPrimVar, val_nPrimVarGrad,
                                                                          config) {

  FcL     = new su2double [nVar];
  FcR     = new su2double [nVar];
  rhos_i  = new su2double [nSpecies];
  rhos_j  = new su2double [nSpecies];
  u_i     = new su2double [nDim];
  u_j     = new su2double [nDim];

  Flux   = new su2double[nVar];

}

CUpwAUSMPWplus_NEMO::~CUpwAUSMPWplus_NEMO(void) {

  delete [] FcL;
  delete [] FcR;
  delete [] rhos_i;
  delete [] rhos_j;
  delete [] u_i;
  delete [] u_j;
  delete [] Flux;
}

CNumerics::ResidualType<> CUpwAUSMPWplus_NEMO::ComputeResidual(const CConfig *config) {

  // NOTE: OSCILLATOR DAMPER "f" NOT IMPLEMENTED!!!

  unsigned short iDim, iVar, iSpecies;
  su2double rho_i, rho_j, rhoEve_i, rhoEve_j;
  su2double aij, atl, gtl_i, gtl_j, sqVi, sqVj, Hnorm;
  su2double w, fL, fR, alpha;
  su2double mL, mR, mLP, mRM, mF, mbLP, mbRM, pLP, pRM, ps;
  su2double gam;

  alpha = 3.0/16.0;

  /*---- Initialize the residual vector ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = 0.0;

  /*--- Calculate geometric quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, c]
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
  rho_i = V_i[RHO_INDEX]; rho_j = V_j[RHO_INDEX];

  rhoCvtr_i = V_i[RHOCVTR_INDEX]; rhoCvtr_j = V_j[RHOCVTR_INDEX];
  rhoCvve_i = V_i[RHOCVVE_INDEX]; rhoCvve_j = V_j[RHOCVVE_INDEX];

  rhoEve_i = 0.0; rhoEve_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhoEve_i += (V_i[RHOS_INDEX+iSpecies]*eve_i[iSpecies]);
    rhoEve_j += (V_j[RHOS_INDEX+iSpecies]*eve_j[iSpecies]);
  }

  /*--- Projected velocities ---*/
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }
  sqVi = 0.0;   sqVj = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqVi += (u_i[iDim]-ProjVel_i*UnitNormal[iDim]) *
            (u_i[iDim]-ProjVel_i*UnitNormal[iDim]);
    sqVj += (u_j[iDim]-ProjVel_j*UnitNormal[iDim]) *
            (u_j[iDim]-ProjVel_j*UnitNormal[iDim]);
  }

  /*--- Calculate interface numerical gammas and speed of sound ---*/
  Hnorm = 0.5*(h_i-0.5*sqVi + h_j-0.5*sqVj);
  gtl_i = Gamma_i;
  gtl_j = Gamma_j;
  gam   = 0.5*(gtl_i+gtl_j);
  if (fabs(rho_i-rho_j)/(0.5*(rho_i+rho_j)) < 1E-3)
    atl = sqrt(2.0*Hnorm*(gam-1.0)/(gam+1.0));
  else {
    atl = sqrt(2.0*Hnorm * (((gtl_i-1.0)/(gtl_i*rho_i) - (gtl_j-1.0)/(gtl_j*rho_j))/
                            ((gtl_j+1.0)/(gtl_j*rho_i) - (gtl_i+1.0)/(gtl_i*rho_j))));
  }

  if (0.5*(ProjVel_i+ProjVel_j) >= 0.0) aij = atl*atl/max(fabs(ProjVel_i),atl);
  else                                  aij = atl*atl/max(fabs(ProjVel_j),atl);

  /*--- Calculate L/R Mach & Pressure functions ---*/
  mL  = ProjVel_i/aij;
  mR  = ProjVel_j/aij;
  if (fabs(mL) <= 1.0) {
    mLP = 0.25*(mL+1.0)*(mL+1.0);
    pLP = P_i*(0.25*(mL+1.0)*(mL+1.0)*(2.0-mL)+alpha*mL*(mL*mL-1.0)*(mL*mL-1.0));
  } else {
    mLP = 0.5*(mL+fabs(mL));
    pLP = P_i*0.5*(mL+fabs(mL))/mL;
  }
  if (fabs(mR) <= 1.0) {
    mRM = -0.25*(mR-1.0)*(mR-1.0);
    pRM = P_j*(0.25*(mR-1.0)*(mR-1.0)*(2.0+mR)-alpha*mR*(mR*mR-1.0)*(mR*mR-1.0));
  } else {
    mRM = 0.5*(mR-fabs(mR));
    pRM = 0.5*P_j*(mR-fabs(mR))/mR;
  }

  /*--- Calculate supporting w & f functions ---*/
  w  = 1.0 - pow(min(P_i/P_j, P_j/P_i), 3.0);
  ps = pLP + pRM;

  // simplified f function (Literature requires information from cells
  // above and below  (TODO)
  if (fabs(mL) < 1.0) fL = P_i/ps - 1.0;
  else fL = 0.0;
  if (fabs(mR) < 1.0) fR = P_j/ps - 1.0;
  else fR = 0.0;

  /*--- Calculate modified M functions ---*/
  mF = mLP + mRM;
  if (mF >= 0.0) {
    mbLP = mLP + mRM*((1.0-w)*(1.0+fR) - fL);
    mbRM = mRM*w*(1.0+fR);
  } else {
    mbLP = mLP*w*(1+fL);
    mbRM = mRM + mLP*((1.0-w)*(1.0+fL) + fL -fR);
  }

  /*--- Assign left & right convective vectors ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    FcL[iSpecies] = rhos_i[iSpecies];
    FcR[iSpecies] = rhos_j[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    FcL[nSpecies+iDim] = rho_i*u_i[iDim];
    FcR[nSpecies+iDim] = rho_j*u_j[iDim];
  }
  FcL[nSpecies+nDim]   = rho_i*h_i;
  FcR[nSpecies+nDim]   = rho_j*h_j;
  FcL[nSpecies+nDim+1] = rhoEve_i;
  FcR[nSpecies+nDim+1] = rhoEve_j;

  /*--- Calculate the numerical flux ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = (mbLP*aij*FcL[iVar] + mbRM*aij*FcR[iVar])*Area;
  for (iDim = 0; iDim < nDim; iDim++)
    Flux[nSpecies+iDim] += (pLP*UnitNormal[iDim] + pRM*UnitNormal[iDim])*Area;

//  if (implicit) //{
//
//    /*--- Initialize the Jacobians ---*/
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_i[iVar][jVar] = 0.0;
//        val_Jacobian_j[iVar][jVar] = 0.0;
//      }
//    }
//
//    /*--- Derivatives of the interface speed of sound, aij ---*/
//    // Derivatives of Hnorm
//    //fact = 0.5*sqrt(2*(gam-1.0)/((gam+1.0)*Hnorm));
//    //for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
//    //  dHnL[iSpecies] = 0.5*(dPdU_i[iSpecies] /*+ sqVi/rho_i*/);
//    //  dHnR[iSpecies] = 0.5*(dPdU_j[iSpecies] /*+ sqVj/rho_j*/);
//    //}
//    //for (iDim = 0; iDim < nDim; iDim++) {
//    //  dV2L = 0.0;
//    //  dV2R = 0.0;
//    //  for (jDim = 0; jDim < nDim; jDim++) {
//    //    dV2L += 2.0/rho_i*(u_i[jDim]-ProjVel_i*UnitNormal[jDim]*(-UnitNormal[iDim]*UnitNormal[jDim]));
//    //    dV2R += 2.0/rho_j*(u_j[jDim]-ProjVel_j*UnitNormal[jDim]*(-UnitNormal[iDim]*UnitNormal[jDim]));
//    //  }
//    //  dV2L += 2.0/rho_i*(u_i[iDim]-ProjVel_i*UnitNormal[iDim] - sqVi);
//    //  dV2R += 2.0/rho_j*(u_j[iDim]-ProjVel_j*UnitNormal[iDim] - sqVj);
//    //  dHnL[nSpecies+iDim] = 0.5*(dPdU_i[nSpecies+iDim] /*- 0.5*(dV2L)*/);
//    //  dHnR[nSpecies+iDim] = 0.5*(dPdU_j[nSpecies+iDim] /*- 0.5*(dV2R)*/);
//    //}
//    //dHnL[nSpecies+nDim]   = 0.5*(1.0+dPdU_i[nSpecies+nDim]);
//    //dHnR[nSpecies+nDim]   = 0.5*(1.0+dPdU_j[nSpecies+nDim]);
//    //dHnL[nSpecies+nDim+1] = 0.5*dPdU_i[nSpecies+nDim+1];
//    //dHnR[nSpecies+nDim+1] = 0.5*dPdU_j[nSpecies+nDim+1];
//
//    //    //////////////////
//    //    //debug:
//    //    cout << "sqVi before: " << sqVi << endl;
//    //    //check sqV routine w/ conserved:
//    //    double rVi, delta;
//    //    rVi = 0.0;
//    //    for (iDim = 0; iDim < nDim; iDim++) {
//    //      rVi += rho_i*u_i[iDim]*UnitNormal[iDim];
//    //    }
//    //    sqVi = 0.0;
//    //    for (iDim = 0; iDim < nDim; iDim++) {
//    //      sqVi += (rho_i*u_i[iDim]-rVi*UnitNormal[iDim])
//    //            * (rho_i*u_i[iDim]-rVi*UnitNormal[iDim])/(rho_i*rho_i);
//    //    }
//    //    cout << "sqVi after: " << sqVi << endl;
//    //
//    //      //perturb:
//    //    delta = V_i[0];
//    //    rho_i = V_i[0]+V_i[1]+delta;
//    //    rVi = 0.0;
//    //    for (iDim = 0; iDim < nDim; iDim++) {
//    //      rVi += rho_i*u_i[iDim]*UnitNormal[iDim];
//    //    }
//    //    sqVj = 0.0;
//    //    for (iDim = 0; iDim < nDim; iDim++) {
//    //      sqVj += (rho_i*u_i[iDim]-rVi*UnitNormal[iDim])
//    //            * (rho_i*u_i[iDim]-rVi*UnitNormal[iDim])/(rho_i*rho_i);
//    //    }
//    //    cout << "FD: " << (sqVj-sqVi)/delta << endl;
//    //    cout << "analytic: " << -2*sqVi/(rho_i-delta) << endl;
//    //    cout << "V0: " << V_i[0] << endl;
//    //    cout << "V1: " << V_i[1] << endl;
//    //    cout << "rho_i: " << rho_i << endl;
//    //    cout << "delta: " << delta << endl;
//    //    cout << "diff: " << sqVj-sqVi << endl;
//    //    cin.get();
//
//
//
//
//    // Derivatives of aij
//    //if (0.5*(ProjVel_i+ProjVel_j) >= 0.0) {
//    //  if (atl >= fabs(ProjVel_i)) {
//    //    for (iVar = 0; iVar < nVar; iVar++) {
//    //      daL[iVar] = fact*dHnL[iVar];
//    //      daR[iVar] = fact*dHnR[iVar];
//    //    }
//    //  } else {
//    //    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    //      daL[iSpecies] = atl*atl/(rho_i*fabs(ProjVel_i))
//    //                    + 2*atl/fabs(ProjVel_i)*fact*dHnL[iSpecies];
//    //      daR[iSpecies] = 2*atl/fabs(ProjVel_i)*fact*dHnR[iSpecies];
//    //    }
//    //    for (iDim = 0; iDim < nDim; iDim++) {
//    //      daL[nSpecies+iDim] = -UnitNormal[iDim]*atl*atl/(fabs(ProjVel_i)*ProjVel_i)
//    //                          + 2*atl/fabs(ProjVel_i)*fact*dHnL[nSpecies+iDim];
//    //       daR[nSpecies+iDim] = 2*atl/fabs(ProjVel_i)*fact*dHnR[nSpecies+iDim];
//    //    }
//    //    daL[nSpecies+nDim]   = 2*atl/fabs(ProjVel_i)*fact*dHnL[nSpecies+nDim];
//    //    daR[nSpecies+nDim]   = 2*atl/fabs(ProjVel_i)*fact*dHnR[nSpecies+nDim];
//    //    daL[nSpecies+nDim+1] = 2*atl/fabs(ProjVel_i)*fact*dHnL[nSpecies+nDim+1];
//    //    daR[nSpecies+nDim+1] = 2*atl/fabs(ProjVel_i)*fact*dHnR[nSpecies+nDim+1];
//    //  }
//    //} else {
//    //  if (atl >= fabs(ProjVel_j)) {
//    //    for (iVar = 0; iVar < nVar; iVar++) {
//    //      daL[iVar] = fact*dHnL[iVar];
//    //      daR[iVar] = fact*dHnR[iVar];
//    //    }
//    //  } else {
//    //    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    //      daR[iSpecies] = atl*atl/(rho_j*fabs(ProjVel_j))
//    //                    + 2*atl/fabs(ProjVel_j)*fact*dHnR[iSpecies];
//    //     daL[iSpecies] = 2*atl/fabs(ProjVel_j)*fact*dHnL[iSpecies];
//    //    }
//    //    for (iDim = 0; iDim < nDim; iDim++) {
//    //      daR[nSpecies+iDim] = -UnitNormal[iDim]*atl*atl/(fabs(ProjVel_j)*ProjVel_j)
//    //                         + 2*atl/fabs(ProjVel_j)*fact*dHnR[nSpecies+iDim];
//    //      daL[nSpecies+iDim] = 2*atl/fabs(ProjVel_j)*fact*dHnL[nSpecies+iDim];
//    //    }
//    //    daR[nSpecies+nDim]   = 2*atl/fabs(ProjVel_j)*fact*dHnR[nSpecies+nDim];
//    //    daL[nSpecies+nDim]   = 2*atl/fabs(ProjVel_j)*fact*dHnL[nSpecies+nDim];
//    //    daR[nSpecies+nDim+1] = 2*atl/fabs(ProjVel_j)*fact*dHnR[nSpecies+nDim+1];
//    //    daL[nSpecies+nDim+1] = 2*atl/fabs(ProjVel_j)*fact*dHnL[nSpecies+nDim+1];
//    //  }
//    // }
//
//    //    cout << "atl: " << atl << endl;
//    //    cout << "ProjVel_i: " << ProjVel_i << endl;
//    //    cout << "term1: " << atl*atl/(rho_i*fabs(ProjVel_i)) << endl;
//    //    cout << "term2: " << endl;
//    //    for (iVar = 0; iVar < nVar; iVar++)
//    //      cout << 2*atl/fabs(ProjVel_i)*fact*dHnL[iVar] << endl;
//    //    cout << "area: " << Area << endl;
//    //    cout << "daL: " << endl;
//    //    for (iVar = 0; iVar < nVar; iVar++) {
//    //      cout << daL[iVar] << endl;
//    //    }
//    //    cin.get();
//
//    /*--- Derivatives of advection speed, mL & mR ---*/
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      dmLdL[iSpecies] = -ProjVel_i/(rho_i*aij) - ProjVel_i/(aij*aij)*daL[iSpecies];
//      dmRdR[iSpecies] = -ProjVel_j/(rho_j*aij) - ProjVel_j/(aij*aij)*daR[iSpecies];
//    }
//    for (iDim = 0; iDim < nDim; iDim++) {
//      dmLdL[nSpecies+iDim] = UnitNormal[iDim]/(rho_i*aij) - ProjVel_i/(aij*aij)*daL[nSpecies+iDim];
//      dmRdR[nSpecies+iDim] = UnitNormal[iDim]/(rho_j*aij) - ProjVel_j/(aij*aij)*daR[nSpecies+iDim];
//    }
//    dmLdL[nSpecies+nDim]   = -ProjVel_i/(aij*aij)*daL[nSpecies+nDim];
//    dmRdR[nSpecies+nDim]   = -ProjVel_j/(aij*aij)*daR[nSpecies+nDim];
//    dmLdL[nSpecies+nDim+1] = -ProjVel_i/(aij*aij)*daL[nSpecies+nDim+1];
//    dmRdR[nSpecies+nDim+1] = -ProjVel_j/(aij*aij)*daR[nSpecies+nDim+1];
//    for (iVar = 0; iVar < nVar; iVar++) {
//      dmLdR[iVar] = -ProjVel_i/(aij*aij)*daR[iVar];
//      dmRdL[iVar] = -ProjVel_j/(aij*aij)*daL[iVar];
//    }
//
//    /*--- Derivatives of numerical advection, mLP & mRM ---*/
//    if (fabs(mL) <= 1.0) {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dmLPdL[iVar] = 0.5*(mL+1)*dmLdL[iVar];
//        dmLPdR[iVar] = 0.5*(mL+1)*dmLdR[iVar];
//      }
//    } else {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dmLPdL[iVar] = 0.5*(dmLdL[iVar] + mL/fabs(mL)*dmLdL[iVar]);
//        dmLPdR[iVar] = 0.5*(dmLdR[iVar] + mL/fabs(mL)*dmLdR[iVar]);
//      }
//    }
//    if (fabs(mR) <= 1.0) {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dmRMdR[iVar] = -0.5*(mR-1)*dmRdR[iVar];
//        dmRMdL[iVar] = -0.5*(mR-1)*dmRdL[iVar];
//      }
//    } else {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dmRMdR[iVar] = 0.5*(dmRdR[iVar] - mR/fabs(mR)*dmRdR[iVar]);
//        dmRMdL[iVar] = 0.5*(dmRdL[iVar] - mR/fabs(mR)*dmRdL[iVar]);
//      }
//    }
//
//    /*--- Derivatives of numerical advection, mbLP & mbRM ---*/
//    if (mF >= 0) {
//      dmbLPdL[iVar] = dmLPdL[iVar] + dmRMdL[iVar]*((1-w)*(1+fR)-fL);
//      dmbLPdR[iVar] = dmLPdR[iVar] + dmRMdR[iVar]*((1-w)*(1+fR)-fL);
//      dmbRMdR[iVar] = dmRMdR[iVar]*w*(1+fR);
//      dmbRMdL[iVar] = dmRMdL[iVar]*w*(1+fR);
//    } else {
//      dmbLPdL[iVar] = dmLPdL[iVar]*w*(1+fL);
//      dmbLPdR[iVar] = dmLPdR[iVar]*w*(1+fL);
//      dmbRMdR[iVar] = dmRMdR[iVar] + dmLPdR[iVar]*((1-w)*(1+fL)+fL-fR);
//      dmbRMdL[iVar] = dmRMdL[iVar] + dmLPdL[iVar]*((1-w)*(1+fL)+fL-fR);
//    }
//
//    /*--- Derivatives of pressure function ---*/
//    if (fabs(mL) <= 1.0) {
//      fact = 0.5*(mL+1)*(2-mL) - 0.25*(mL+1)*(mL+1)
//          + alpha*(mL*mL-1)*(mL*mL-1) + 4*alpha*mL*mL*(mL*mL-1);
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dpLPdL[iVar] = dPdU_i[iVar]*pLP/P_i + P_i*fact*dmLdL[iVar];
//        dpLPdR[iVar] = P_i*fact*dmLdR[iVar];
//      }
//    } else {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dpLPdL[iVar] = dPdU_i[iVar] * 0.5*(mL+fabs(mL))/mL;
//        dpLPdR[iVar] = 0.0;
//      }
//    }
//    if (fabs(mR) <= 1.0) {
//      fact = 0.5*(mR-1)*(2+mR) + 0.25*(mR-1)*(mR-1)
//          - alpha*(mR*mR-1)*(mR*mR-1) - 4*alpha*mR*mR*(mR*mR-1);
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dpRMdR[iVar] = dPdU_j[iVar]*pRM/P_j + P_j*fact*dmRdR[iVar];
//        dpRMdL[iVar] = P_j*fact*dmRdL[iVar];
//      }
//    } else {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dpRMdR[iVar] = dPdU_j[iVar] * 0.5*(mR+fabs(mR))/mR;
//        dpRMdL[iVar] = 0.0;
//      }
//    }
//
//    /*--- L Jacobian ---*/
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_i[iVar][jVar] += (dmbLPdL[jVar]*FcL[iVar] + dmbRMdL[jVar]*FcR[iVar])*aij*Area;
//        val_Jacobian_i[iVar][jVar] += (mbLP*FcL[iVar] + mbRM*FcR[iVar])*daL[jVar]*Area;
//      }
//      val_Jacobian_i[iVar][iVar] += mbLP*aij*Area;
//      val_Jacobian_i[nSpecies+nDim][iVar] += mbLP*aij*dPdU_i[iVar]*Area;
//
//      // pressure terms
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_Jacobian_i[nSpecies+iDim][iVar] += dpLPdL[iVar]*UnitNormal[iDim]*Area;
//        val_Jacobian_i[nSpecies+iDim][iVar] += dpRMdL[iVar]*UnitNormal[iDim]*Area;
//      }
//    }
//    /*--- R Jacobian ---*/
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_j[iVar][jVar] += (dmbLPdR[jVar]*FcL[iVar] + dmbRMdR[jVar]*FcR[iVar])*aij*Area;
//        val_Jacobian_j[iVar][jVar] += (mbLP*FcL[iVar] + mbRM*FcR[iVar])*daR[jVar]*Area;
//      }
//      val_Jacobian_j[iVar][iVar] += mbRM*aij*Area;
//      val_Jacobian_j[nSpecies+nDim][iVar] += mbRM*aij*dPdU_j[iVar]*Area;
//
//      // pressure terms
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_Jacobian_j[nSpecies+iDim][iVar] += dpLPdR[iVar]*UnitNormal[iDim]*Area;
//        val_Jacobian_j[nSpecies+iDim][iVar] += dpRMdR[iVar]*UnitNormal[iDim]*Area;
//      }
//    }
//  }

  return ResidualType<>(Flux, nullptr, nullptr);
}
