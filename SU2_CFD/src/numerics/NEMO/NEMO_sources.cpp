/*!
 * \file NEMO_sources.cpp
 * \brief Implementation of numerics classes for integration
 *        of source terms in fluid flow NEMO problems.
 * \author C. Garbacz, W. Maier, S. Copeland.
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/numerics/NEMO/NEMO_sources.hpp"

CSource_NEMO::CSource_NEMO(unsigned short val_nDim,
                           unsigned short val_nVar,
                           unsigned short val_nPrimVar,
                           unsigned short val_nPrimVarGrad,
                           CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  unsigned short iSpecies;

  /*--- Allocate arrays ---*/
  alphak = new int[nSpecies];
  betak  = new int[nSpecies];
  Y      = new su2double[nSpecies];
  dkf    = new su2double[nVar];
  dkb    = new su2double[nVar];
  dRfok  = new su2double[nVar];
  dRbok  = new su2double[nVar];

  ws.resize(nSpecies,0.0);

  dYdr = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dYdr[iSpecies] = new su2double[nSpecies];
  }

  residual = new su2double[nVar];
}

CSource_NEMO::~CSource_NEMO(void) {
  unsigned short iSpecies;

  /*--- Deallocate arrays ---*/

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] dYdr[iSpecies];
  delete [] dYdr;

  delete [] Y;
  delete [] alphak;
  delete [] betak;
  delete [] dkf;
  delete [] dkb;
  delete [] dRfok;
  delete [] dRbok;

  delete [] residual;

}

CNumerics::ResidualType<> CSource_NEMO::ComputeChemistry(const CConfig *config) {

  /*--- Nonequilibrium chemistry ---*/
  unsigned short iSpecies, iVar;
  su2double T, Tve;
  vector<su2double> rhos;

  rhos.resize(nSpecies,0.0);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }

  /*--- Rename for convenience ---*/
  T       = V_i[T_INDEX];
  Tve     = V_i[TVE_INDEX];
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies]=V_i[RHOS_INDEX+iSpecies];

  /*--- Set mixture state ---*/
  fluidmodel->SetTDStateRhosTTv(rhos, T, Tve);

  ws = fluidmodel->ComputeNetProductionRates();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    residual[iSpecies] = ws[iSpecies] *Volume;

  //if (implicit) {
  //  for (iVar = 0; iVar < nVar; iVar++)
  //    for (jVar = 0; jVar < nVar; jVar++)
  //      val_Jacobian_i[iVar][jVar] = 0.0;
  //}
//  if (implicit) {
//    su2double dThf, dThb;
//
//      /*--- Initializing derivative variables ---*/
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dkf[iVar] = 0.0;
//        dkb[iVar] = 0.0;
//        dRfok[iVar] = 0.0;
//        dRbok[iVar] = 0.0;
//      }
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        alphak[iSpecies] = 0;
//        betak[iSpecies]  = 0;
//      }
//
//      /*--- Derivative of modified temperature wrt Trxnf ---*/
//      dThf = 0.5 * (1.0 + (Trxnf-T_min)/sqrt((Trxnf-T_min)*(Trxnf-T_min)
//                                             + epsilon*epsilon          ));
//      dThb = 0.5 * (1.0 + (Trxnb-T_min)/sqrt((Trxnb-T_min)*(Trxnb-T_min)
//                                             + epsilon*epsilon          ));
//
//      /*--- Fwd rate coefficient derivatives ---*/
//      coeff = kf * (eta/Thf+theta/(Thf*Thf)) * dThf;
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dkf[iVar] = coeff * (  af*Trxnf/T*dTdU_i[iVar]
//                               + bf*Trxnf/Tve*dTvedU_i[iVar] );
//      }
//
//      /*--- Bkwd rate coefficient derivatives ---*/
//      coeff = kb * (eta/Thb+theta/(Thb*Thb)) * dThb;
//      for (iVar = 0; iVar < nVar; iVar++) {
//        dkb[iVar] = coeff*(  ab*Trxnb/T*dTdU_i[iVar]
//                             + bb*Trxnb/Tve*dTvedU_i[iVar])
//            - kb*((A[0]*Thb/1E4 - A[2] - A[3]*1E4/Thb
//            - 2*A[4]*(1E4/Thb)*(1E4/Thb))/Thb) * dThb * (  ab*Trxnb/T*dTdU_i[iVar]
//                                                           + bb*Trxnb/Tve*dTvedU_i[iVar]);
//      }
//
//      /*--- Rxn rate derivatives ---*/
//      for (ii = 0; ii < 3; ii++) {
//
//        /*--- Products ---*/
//        iSpecies = RxnMap(iReaction,1,ii);
//        if (iSpecies != nSpecies)
//          betak[iSpecies]++;
//
//        /*--- Reactants ---*/
//        iSpecies = RxnMap(iReaction,0,ii);
//        if (iSpecies != nSpecies)
//          alphak[iSpecies]++;
//      }
//
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//
//        // Fwd
//        dRfok[iSpecies] =  0.001*alphak[iSpecies]/Ms[iSpecies]
//            * pow(0.001*U_i[iSpecies]/Ms[iSpecies],
//                  max(0, alphak[iSpecies]-1)      );
//        for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
//          if (jSpecies != iSpecies)
//            dRfok[iSpecies] *= pow(0.001*U_i[jSpecies]/Ms[jSpecies],
//                                   alphak[jSpecies]                );
//        dRfok[iSpecies] *= 1000.0;
//
//        // Bkw
//        dRbok[iSpecies] =  0.001*betak[iSpecies]/Ms[iSpecies]
//            * pow(0.001*U_i[iSpecies]/Ms[iSpecies],
//                  max(0, betak[iSpecies]-1)       );
//        for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
//          if (jSpecies != iSpecies)
//            dRbok[iSpecies] *= pow(0.001*U_i[jSpecies]/Ms[jSpecies],
//                                   betak[jSpecies]                 );
//        dRbok[iSpecies] *= 1000.0;
//      }
//
//      nEve = nSpecies+nDim+1;
//      for (ii = 0; ii < 3; ii++) {
//
//        /*--- Products ---*/
//        iSpecies = RxnMap(iReaction,1,ii);
//        if (iSpecies != nSpecies) {
//          for (iVar = 0; iVar < nVar; iVar++) {
//            val_Jacobian_i[iSpecies][iVar] +=
//                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
//                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) * Volume;
//            val_Jacobian_i[nEve][iVar] +=
//                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
//                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar])
//                * eve_i[iSpecies] * Volume;
//          }
//
//          for (jVar = 0; jVar < nVar; jVar++) {
//            val_Jacobian_i[nEve][jVar] += Ms[iSpecies] * (fwdRxn-bkwRxn)
//                * Cvve_i[iSpecies] * dTvedU_i[jVar] * Volume;
//          }
//        }
//
//        /*--- Reactants ---*/
//        iSpecies = RxnMap(iReaction,0,ii);
//        if (iSpecies != nSpecies) {
//          for (iVar = 0; iVar < nVar; iVar++) {
//            val_Jacobian_i[iSpecies][iVar] -=
//                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
//                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) * Volume;
//            val_Jacobian_i[nEve][iVar] -=
//                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
//                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar])
//                * eve_i[iSpecies] * Volume;
//
//          }
//
//          for (jVar = 0; jVar < nVar; jVar++) {
//            val_Jacobian_i[nEve][jVar] -= Ms[iSpecies] * (fwdRxn-bkwRxn)
//                * Cvve_i[iSpecies] * dTvedU_i[jVar] * Volume;
//          }
//        } // != nSpecies
//      } // ii
//    } // implicit

  return ResidualType<>(residual, nullptr, nullptr);

}

CNumerics::ResidualType<> CSource_NEMO::ComputeVibRelaxation(const CConfig *config) {

  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section
  unsigned short iSpecies, iVar;
  su2double  T, Tve;
  su2double res_min = -1E6;
  su2double res_max = 1E6;
  vector<su2double> rhos;

  rhos.resize(nSpecies,0.0);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }

  /*--- Rename for convenience ---*/
  T       = V_i[T_INDEX];
  Tve     = V_i[TVE_INDEX];
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies]=V_i[RHOS_INDEX+iSpecies];

  fluidmodel->SetTDStateRhosTTv(rhos, T, Tve);

  residual[nSpecies+nDim+1] = fluidmodel->ComputeEveSourceTerm() * Volume;

  //  if (implicit) {
//    for (iVar = 0; iVar < nVar; iVar++)
//      for (jVar = 0; jVar < nVar; jVar++)
//        val_Jacobian_i[iVar][jVar] = 0.0;
//  }

//  if (implicit) {
//
//    fluidmodel->SetTve(T);
//    Cvvsst = fluidmodel->ComputeSpeciesCvVibEle();
//
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//
//      for (iVar = 0; iVar < nVar; iVar++) {
//        val_Jacobian_i[nSpecies+nDim+1][iVar] += rhos[iSpecies]/taus[iSpecies]*(Cvvsst[iSpecies]*dTdU_i[iVar] -
//                                                          Cvve_i[iSpecies]*dTvedU_i[iVar])*Volume;
//      }
//    }
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      val_Jacobian_i[nSpecies+nDim+1][iSpecies] += (estar[iSpecies]-eve_i[iSpecies])/taus[iSpecies]*Volume;
//  }

  if(config->GetVTTransferResidualLimiting()){
    if(residual[nSpecies+nDim+1]>res_max) residual[nSpecies+nDim+1]=res_max;
    if(residual[nSpecies+nDim+1]<res_min) residual[nSpecies+nDim+1]=res_min;
  }

  return ResidualType<>(residual, nullptr, nullptr);
}

CNumerics::ResidualType<> CSource_NEMO::ComputeAxisymmetric(const CConfig *config) {

  unsigned short iDim, iSpecies, iVar;
  su2double rho, rhov, vel2, H, yinv, T, Tve, Ru, RuSI;
  su2double ktr, kve;

  /*--- Rename for convenience ---*/
  auto Ds = Diffusion_Coeff_i;
  ktr = Thermal_Conductivity_i;
  kve = Thermal_Conductivity_ve_i;
  rho = V_i[RHO_INDEX];
  T   = V_i[T_INDEX];
  Tve  = V_i[TVE_INDEX];
  auto GV  = PrimVar_Grad_i;
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;

  auto& Ms = fluidmodel->GetSpeciesMolarMass();

  bool viscous = config->GetViscous();
  bool rans = (config->GetKind_Turb_Model() != NONE);
  hs = fluidmodel->ComputeSpeciesEnthalpy(T, Tve, eve_i);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }

  /*--- Calculate inverse of y coordinate ---*/
  if (Coord_i[1]!= 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;

  /*--- Rename for convenience ---*/
  rho    = V_i[RHO_INDEX];
  rhov   = U_i[nSpecies+1];
  H      = V_i[H_INDEX];
  vel2   = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    vel2 += V_i[VEL_INDEX+iDim]*V_i[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Y[iSpecies] = V_i[RHOS_INDEX+iSpecies] / rho;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    residual[iSpecies] = yinv*rhov*Y[iSpecies]*Volume;
  residual[nSpecies]   = yinv*rhov*U_i[nSpecies]/rho*Volume;
  residual[nSpecies+1] = yinv*rhov*U_i[nSpecies+1]/rho*Volume;
  residual[nSpecies+2] = yinv*rhov*H*Volume;
  residual[nSpecies+3] = yinv*rhov*U_i[nSpecies+nDim+1]/rho*Volume;

  if (viscous) {
    if (!rans){ turb_ke_i = 0.0; }

    su2double Vector = 0.0;
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
      Vector += rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][1];

    su2double sumJhs_y  = 0.0;
    su2double sumJeve_y = 0.0;
    su2double Mass      = 0.0;

    for (iSpecies=0; iSpecies<nSpecies; iSpecies++)
      Mass += V_i[iSpecies]*Ms[iSpecies];

    su2double heat_capacity_cp_i   = V_i[RHOCVTR_INDEX]+Ru/Mass;
    su2double total_viscosity_i    = Laminar_Viscosity_i + Eddy_Viscosity_i;
    su2double total_conductivity_i = ktr + kve + heat_capacity_cp_i*Eddy_Viscosity_i/Prandtl_Turb;
    su2double u                    = V_i[VEL_INDEX];
    su2double v                    = V_i[VEL_INDEX+1];
    su2double qy_ve                = kve*GV[TVE_INDEX][1];

    /*--- Enthalpy and vib-el energy transport due to y-direction diffusion---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      sumJhs_y  += (rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][1] - V_i[RHOS_INDEX+iSpecies]*Vector) * hs[iSpecies];
      sumJeve_y += (rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][1] - V_i[RHOS_INDEX+iSpecies]*Vector) * eve_i[iSpecies];
    }

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      residual[iSpecies] -= 0.0;
    residual[nSpecies] -= Volume*(yinv*total_viscosity_i*(PrimVar_Grad_i[nSpecies+2][1]+PrimVar_Grad_i[nSpecies+3][0])
                                                                    -TWO3*AuxVar_Grad_i[0][0]);
    residual[nSpecies+1] -= Volume*(yinv*total_viscosity_i*2*(PrimVar_Grad_i[nSpecies+3][1]-v*yinv)
                                                                    -TWO3*AuxVar_Grad_i[0][1]);
    residual[nSpecies+2] -= Volume*(yinv*(sumJhs_y + total_viscosity_i*(u*(PrimVar_Grad_i[nSpecies+3][0]+PrimVar_Grad_i[nSpecies+2][1])
                                                     +v*TWO3*(2*PrimVar_Grad_i[nSpecies+2][1]-PrimVar_Grad_i[nSpecies+2][0]
                                                     -v*yinv+rho*turb_ke_i))
                                                     -total_conductivity_i*PrimVar_Grad_i[nSpecies][1])
                                                     -TWO3*(AuxVar_Grad_i[1][1]+AuxVar_Grad_i[2][0]));
    residual[nSpecies+3] -= Volume*(yinv*(sumJeve_y -qy_ve));
  }

//  if (implicit) {
//
//    /*--- Initialize ---*/
//    for (iVar = 0; iVar < nVar; iVar++)
//      for (jVar = 0; jVar < nVar; jVar++)
//        val_Jacobian[iVar][jVar] = 0.0;
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
//        dYdr[iSpecies][jSpecies] = 0.0;
//
//    /*--- Calculate additional quantities ---*/
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//        dYdr[iSpecies][jSpecies] += -1/rho*Ys[iSpecies];
//      }
//      dYdr[iSpecies][iSpecies] += 1/rho;
//    }
//
//    /*--- Populate Jacobian ---*/
//
//    // Species density
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//        val_Jacobian[iSpecies][jSpecies] = dYdr[iSpecies][jSpecies]*rhov;
//      }
//      val_Jacobian[iSpecies][nSpecies+1] = Y[iSpecies];
//    }
//
//    // X-momentum
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      val_Jacobian[nSpecies][iSpecies] = -rhou*rhov/(rho*rho);
//    val_Jacobian[nSpecies][nSpecies] = rhov/rho;
//    val_Jacobian[nSpecies][nSpecies+1] = rhou/rho;
//
//    // Y-momentum
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      val_Jacobian[nSpecies+1][iSpecies] = -rhov*rhov/(rho*rho);
//    val_Jacobian[nSpecies+1][nSpecies+1] = 2*rhov/rho;
//
//    // Energy
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      val_Jacobian[nSpecies+nDim][iSpecies]      = -H*rhov/rho + dPdU_i[iSpecies]*rhov/rho;
//    val_Jacobian[nSpecies+nDim][nSpecies]        = dPdU_i[nSpecies]*rhov/rho;
//    val_Jacobian[nSpecies+nDim][nSpecies+1]      = H + dPdU_i[nSpecies+1]*rhov/rho;
//    val_Jacobian[nSpecies+nDim][nSpecies+nDim]   = (1+dPdU_i[nSpecies+nDim])*rhov/rho;
//    val_Jacobian[nSpecies+nDim][nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1]*rhov/rho;
//
//    // Vib-el energy
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      val_Jacobian[nSpecies+nDim+1][iSpecies] = -rhoEve*rhov/(rho*rho);
//    val_Jacobian[nSpecies+nDim+1][nSpecies+1] = rhoEve/rho;
//    val_Jacobian[nSpecies+nDim+1][nSpecies+nDim+1] = rhov/rho;
//
//    for (iVar = 0; iVar < nVar; iVar++)
//      for (jVar = 0; jVar < nVar; jVar++)
//        val_Jacobian[iVar][jVar] *= yinv*Volume;
//  }

  return ResidualType<>(residual, nullptr, nullptr);
}

