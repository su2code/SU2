/*!
 * \file NEMO_sources.cpp
 * \brief Implementation of numerics classes for integration
 *        of source terms in fluid flow NEMO problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
                           CConfig *config) : CNumerics(val_nDim,
                                                        val_nVar,
                                                        config) {

  unsigned short iVar, iSpecies;

  /*--- Assign booleans from CConfig ---*/
  implicit   = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  ionization = config->GetIonization();

  /*--- Define useful constants ---*/
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();

  /*--- Allocate arrays ---*/
  RxnConstantTable = new su2double*[6];
  for (iVar = 0; iVar < 6; iVar++)
    RxnConstantTable[iVar] = new su2double[5];

  tau_sr = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    tau_sr[iSpecies] = new su2double[nSpecies];

  alphak = new int[nSpecies];
  betak  = new int[nSpecies];
  A      = new su2double[5];
  X      = new su2double[nSpecies];
  Y      = new su2double[nSpecies];
  estar  = new su2double[nSpecies];
  evib   = new su2double[nSpecies];
  Cvvs   = new su2double[nSpecies];
  Cves   = new su2double[nSpecies];
  Cvvsst = new su2double[nSpecies];
  tauP   = new su2double[nSpecies];
  tauMW  = new su2double[nSpecies];
  taus   = new su2double[nSpecies];
  dkf    = new su2double[nVar];
  dkb    = new su2double[nVar];
  dRfok  = new su2double[nVar];
  dRbok  = new su2double[nVar];

  dYdr = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dYdr[iSpecies] = new su2double[nSpecies];
  }

  variable = new CNEMOEulerVariable(1, nDim, nVar, nPrimVar, nPrimVarGrad, config);
}

CSource_NEMO::~CSource_NEMO(void) {
  unsigned short iVar, iSpecies;

  /*--- Deallocate arrays ---*/

  for (iVar = 0; iVar < 6; iVar++)
    delete [] RxnConstantTable[iVar];
  delete [] RxnConstantTable;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] dYdr[iSpecies];
  delete [] dYdr;

  for (iSpecies = 0;iSpecies < nSpecies; iSpecies++)
    delete [] tau_sr[iSpecies];
  delete [] tau_sr;

  delete [] A;
  delete [] X;
  delete [] Y;
  delete [] estar;
  delete [] evib;
  delete [] Cvvs;
  delete [] Cves;
  delete [] Cvvsst;
  delete [] tauP;
  delete [] tauMW;
  delete [] taus;
  delete [] alphak;
  delete [] betak;
  delete [] dkf;
  delete [] dkb;
  delete [] dRfok;
  delete [] dRbok;

}

void CSource_NEMO::GetKeqConstants(su2double *A, unsigned short val_Reaction,
                                   CConfig *config) {
  unsigned short ii, iSpecies, iIndex, tbl_offset, pwr;
  su2double N;
  su2double tmp1, tmp2;

  /*--- Acquire database constants from CConfig ---*/
  const su2double *Ms = config->GetMolar_Mass();
  config->GetChemistryEquilConstants(RxnConstantTable, val_Reaction);

  /*--- Calculate mixture number density ---*/
  N = 0.0;
  for (iSpecies =0 ; iSpecies < nSpecies; iSpecies++) {
    N += V_i[iSpecies]/Ms[iSpecies]*AVOGAD_CONSTANT;
  }

  /*--- Convert number density from 1/m^3 to 1/cm^3 for table look-up ---*/
  N = N*(1E-6);

  /*--- Determine table index based on mixture N ---*/
  tbl_offset = 14;
  pwr        = floor(log10(N));

  /*--- Bound the interpolation to table limit values ---*/
  iIndex = int(pwr) - tbl_offset;
  if (iIndex <= 0) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable[0][ii];
    return;
  } else if (iIndex >= 5) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable[5][ii];
    return;
  }

  /*--- Calculate interpolation denominator terms avoiding pow() ---*/
  tmp1 = 1.0;
  tmp2 = 1.0;
  for (ii = 0; ii < pwr; ii++) {
    tmp1 *= 10.0;
    tmp2 *= 10.0;
  }
  tmp2 *= 10.0;

  /*--- Interpolate ---*/
  for (ii = 0; ii < 5; ii++) {
    A[ii] =  (RxnConstantTable[iIndex+1][ii] - RxnConstantTable[iIndex][ii])
        / (tmp2 - tmp1) * (N - tmp1)
        + RxnConstantTable[iIndex][ii];
  }
  return;
}

void CSource_NEMO::ComputeChemistry(su2double *val_residual,
                                    su2double *val_source,
                                    su2double **val_Jacobian_i,
                                    CConfig *config) {

  /*--- Nonequilibrium chemistry ---*/
  unsigned short iSpecies, jSpecies, ii, iReaction, nReactions, iVar, jVar;
  unsigned short nHeavy, nEl, nEve;
  su2double T_min, epsilon;
  su2double T, Tve, Thf, Thb, Trxnf, Trxnb, Keq, Cf, eta, theta, kf, kb, kfb;
  su2double rho, rhoCvtr, rhoCvve, P;
  su2double fwdRxn, bkwRxn, alpha, RuSI, Ru;
  su2double af, bf, ab, bb, coeff;
  su2double dThf, dThb;

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
  }

  /*--- Define artificial chemistry parameters ---*/
  // Note: These parameters artificially increase the rate-controlling reaction
  //       temperature.  This relaxes some of the stiffness in the chemistry
  //       source term.
  T_min   = 800.0;
  epsilon = 80;

  /*--- Define preferential dissociation coefficient ---*/
  alpha = 0.3;

  /*--- Determine the number of heavy particle species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Rename for convenience ---*/
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;
//  cout << "delete me : "<< V_i[T_INDEX] << endl;
//  cout << "delete me : "<< V_i[TVE_INDEX] << endl;
//  cout << "delete me : "<< V_i[P_INDEX] << endl;
  rho     = V_i[RHO_INDEX];
  P       = V_i[P_INDEX];
  T       = V_i[T_INDEX];
  Tve     = V_i[TVE_INDEX];
  rhoCvtr = V_i[RHOCVTR_INDEX];
  rhoCvve = V_i[RHOCVVE_INDEX];

  /*--- Acquire parameters from the configuration file ---*/
  nReactions = config->GetnReactions();
  const su2double *Ms         = config->GetMolar_Mass();
  const su2double *hf         = config->GetEnthalpy_Formation();
  const su2double *xi         = config->GetRotationModes();
  const su2double *Tref       = config->GetRefTemperature();
  const su2double *Tcf_a      = config->GetRxnTcf_a();
  const su2double *Tcf_b      = config->GetRxnTcf_b();
  const su2double *Tcb_a      = config->GetRxnTcb_a();
  const su2double *Tcb_b      = config->GetRxnTcb_b();
  const auto& RxnMap = config->GetReaction_Map();

  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    /*--- Determine the rate-controlling temperature ---*/
    af = Tcf_a[iReaction];
    bf = Tcf_b[iReaction];
    ab = Tcb_a[iReaction];
    bb = Tcb_b[iReaction];
    Trxnf = pow(T, af)*pow(Tve, bf);
    Trxnb = pow(T, ab)*pow(Tve, bb);

    /*--- Calculate the modified temperature ---*/
    Thf = 0.5 * (Trxnf+T_min + sqrt((Trxnf-T_min)*(Trxnf-T_min)+epsilon*epsilon));
    Thb = 0.5 * (Trxnb+T_min + sqrt((Trxnb-T_min)*(Trxnb-T_min)+epsilon*epsilon));

    /*--- Get the Keq & Arrhenius coefficients ---*/
    GetKeqConstants(A, iReaction, config);
    Cf    = config->GetArrheniusCoeff(iReaction);
    eta   = config->GetArrheniusEta(iReaction);
    theta = config->GetArrheniusTheta(iReaction);

    /*--- Calculate Keq ---*/
    Keq = exp(  A[0]*(Thb/1E4) + A[1] + A[2]*log(1E4/Thb)
        + A[3]*(1E4/Thb) + A[4]*(1E4/Thb)*(1E4/Thb) );

    /*--- Calculate rate coefficients ---*/
    kf  = Cf * exp(eta*log(Thf)) * exp(-theta/Thf);
    kfb = Cf * exp(eta*log(Thb)) * exp(-theta/Thb);
    kb  = kfb / Keq;

    /*--- Determine production & destruction of each species ---*/
    fwdRxn = 1.0;
    bkwRxn = 1.0;
    for (ii = 0; ii < 3; ii++) {

      /*--- Reactants ---*/
      iSpecies = RxnMap(iReaction,0,ii);
      if ( iSpecies != nSpecies)
        fwdRxn *= 0.001*U_i[iSpecies]/Ms[iSpecies];

      /*--- Products ---*/
      jSpecies = RxnMap(iReaction,1,ii);
      if (jSpecies != nSpecies) {
        bkwRxn *= 0.001*U_i[jSpecies]/Ms[jSpecies];
      }
    }
    fwdRxn = 1000.0 * kf * fwdRxn;
    bkwRxn = 1000.0 * kb * bkwRxn;

    for (ii = 0; ii < 3; ii++) {

      /*--- Products ---*/
      iSpecies = RxnMap(iReaction,1,ii);
      if (iSpecies != nSpecies) {
        val_residual[iSpecies] += Ms[iSpecies] * (fwdRxn-bkwRxn) * Volume;
        val_residual[nSpecies+nDim+1] += Ms[iSpecies] * (fwdRxn-bkwRxn)
            * eve_i[iSpecies] * Volume;
      }

      /*--- Reactants ---*/
      iSpecies = RxnMap(iReaction,0,ii);
      if (iSpecies != nSpecies) {
        val_residual[iSpecies] -= Ms[iSpecies] * (fwdRxn-bkwRxn) * Volume;
        val_residual[nSpecies+nDim+1] -= Ms[iSpecies] * (fwdRxn-bkwRxn)
            * eve_i[iSpecies] * Volume;
      }
    }

    /*---Set source term ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      val_source[iVar] = val_source[iVar]+val_residual[iVar]/Volume;

    if (implicit) {

      /*--- Initializing derivative variables ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        dkf[iVar] = 0.0;
        dkb[iVar] = 0.0;
        dRfok[iVar] = 0.0;
        dRbok[iVar] = 0.0;
      }
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        alphak[iSpecies] = 0;
        betak[iSpecies]  = 0;
      }

      /*--- Derivative of modified temperature wrt Trxnf ---*/
      dThf = 0.5 * (1.0 + (Trxnf-T_min)/sqrt((Trxnf-T_min)*(Trxnf-T_min)
                                             + epsilon*epsilon          ));
      dThb = 0.5 * (1.0 + (Trxnb-T_min)/sqrt((Trxnb-T_min)*(Trxnb-T_min)
                                             + epsilon*epsilon          ));

      /*--- Fwd rate coefficient derivatives ---*/
      coeff = kf * (eta/Thf+theta/(Thf*Thf)) * dThf;
      for (iVar = 0; iVar < nVar; iVar++) {
        dkf[iVar] = coeff * (  af*Trxnf/T*dTdU_i[iVar]
                               + bf*Trxnf/Tve*dTvedU_i[iVar] );
      }

      /*--- Bkwd rate coefficient derivatives ---*/
      coeff = kb * (eta/Thb+theta/(Thb*Thb)) * dThb;
      for (iVar = 0; iVar < nVar; iVar++) {
        dkb[iVar] = coeff*(  ab*Trxnb/T*dTdU_i[iVar]
                             + bb*Trxnb/Tve*dTvedU_i[iVar])
            - kb*((A[0]*Thb/1E4 - A[2] - A[3]*1E4/Thb
            - 2*A[4]*(1E4/Thb)*(1E4/Thb))/Thb) * dThb * (  ab*Trxnb/T*dTdU_i[iVar]
                                                           + bb*Trxnb/Tve*dTvedU_i[iVar]);
      }

      /*--- Rxn rate derivatives ---*/
      for (ii = 0; ii < 3; ii++) {

        /*--- Products ---*/
        iSpecies = RxnMap(iReaction,1,ii);
        if (iSpecies != nSpecies)
          betak[iSpecies]++;

        /*--- Reactants ---*/
        iSpecies = RxnMap(iReaction,0,ii);
        if (iSpecies != nSpecies)
          alphak[iSpecies]++;
      }

      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

        // Fwd
        dRfok[iSpecies] =  0.001*alphak[iSpecies]/Ms[iSpecies]
            * pow(0.001*U_i[iSpecies]/Ms[iSpecies],
                  max(0, alphak[iSpecies]-1)      );
        for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
          if (jSpecies != iSpecies)
            dRfok[iSpecies] *= pow(0.001*U_i[jSpecies]/Ms[jSpecies],
                                   alphak[jSpecies]                );
        dRfok[iSpecies] *= 1000.0;

        // Bkw
        dRbok[iSpecies] =  0.001*betak[iSpecies]/Ms[iSpecies]
            * pow(0.001*U_i[iSpecies]/Ms[iSpecies],
                  max(0, betak[iSpecies]-1)       );
        for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
          if (jSpecies != iSpecies)
            dRbok[iSpecies] *= pow(0.001*U_i[jSpecies]/Ms[jSpecies],
                                   betak[jSpecies]                 );
        dRbok[iSpecies] *= 1000.0;
      }

      nEve = nSpecies+nDim+1;
      for (ii = 0; ii < 3; ii++) {

        /*--- Products ---*/
        iSpecies = RxnMap(iReaction,1,ii);
        if (iSpecies != nSpecies) {
          for (iVar = 0; iVar < nVar; iVar++) {
            val_Jacobian_i[iSpecies][iVar] +=
                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) * Volume;
            val_Jacobian_i[nEve][iVar] +=
                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar])
                * eve_i[iSpecies] * Volume;
          }

          for (jVar = 0; jVar < nVar; jVar++) {
            val_Jacobian_i[nEve][jVar] += Ms[iSpecies] * (fwdRxn-bkwRxn)
                * Cvve_i[iSpecies] * dTvedU_i[jVar] * Volume;
          }
        }

        /*--- Reactants ---*/
        iSpecies = RxnMap(iReaction,0,ii);
        if (iSpecies != nSpecies) {
          for (iVar = 0; iVar < nVar; iVar++) {
            val_Jacobian_i[iSpecies][iVar] -=
                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) * Volume;
            val_Jacobian_i[nEve][iVar] -=
                Ms[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar]
                                 -dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar])
                * eve_i[iSpecies] * Volume;

          }

          for (jVar = 0; jVar < nVar; jVar++) {
            val_Jacobian_i[nEve][jVar] -= Ms[iSpecies] * (fwdRxn-bkwRxn)
                * Cvve_i[iSpecies] * dTvedU_i[jVar] * Volume;
          }
        } // != nSpecies
      } // ii
    } // implicit
  } // iReaction
}

void CSource_NEMO::ComputeVibRelaxation(su2double *val_residual,
                                        su2double *val_source,
                                        su2double **val_Jacobian_i,
                                        CConfig *config) {

  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section
  unsigned short iSpecies, jSpecies, iVar, jVar;
  unsigned short nEv, nHeavy, nEl;
  su2double rhos, P, T, Tve, rhoCvtr, rhoCvve, RuSI, Ru, conc, N;
  su2double Qtv, taunum, taudenom;
  su2double mu, A_sr, B_sr, num, denom;
  su2double Cs, sig_s;

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
  }

  /*--- Determine the number of heavy particle species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Rename for convenience ---*/
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;
  P       = V_i[P_INDEX];
  T       = V_i[T_INDEX];
  Tve     = V_i[TVE_INDEX];
  rhoCvtr = V_i[RHOCVTR_INDEX];
  rhoCvve = V_i[RHOCVVE_INDEX];
  nEv     = nSpecies+nDim+1;

  /*--- Read from CConfig ---*/
  const su2double *Ms        = config->GetMolar_Mass();
  const su2double *thetav    = config->GetCharVibTemp();

  /*--- Calculate mole fractions ---*/
  N    = 0.0;
  conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc += V_i[RHOS_INDEX+iSpecies] / Ms[iSpecies];
    N    += V_i[RHOS_INDEX+iSpecies] / Ms[iSpecies] * AVOGAD_CONSTANT;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    X[iSpecies] = (V_i[RHOS_INDEX+iSpecies] / Ms[iSpecies]) / conc;

  /*--- Loop over species to calculate source term --*/
  Qtv      = 0.0;
  taunum   = 0.0;
  taudenom = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Rename for convenience ---*/
    rhos   = V_i[RHOS_INDEX+iSpecies];

    /*--- Millikan & White relaxation time ---*/
    num   = 0.0;
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      mu     = Ms[iSpecies]*Ms[jSpecies] / (Ms[iSpecies] + Ms[jSpecies]);
      A_sr   = 1.16 * 1E-3 * sqrt(mu) * pow(thetav[iSpecies], 4.0/3.0);
      B_sr   = 0.015 * pow(mu, 0.25);
      tau_sr[iSpecies][jSpecies] = 101325.0/P * exp(A_sr*(pow(T,-1.0/3.0) - B_sr) - 18.42);
      num   += X[jSpecies];
      denom += X[jSpecies] / tau_sr[iSpecies][jSpecies];
    }

    tauMW[iSpecies] = num / denom;
  
    /*--- Park limiting cross section ---*/
    Cs    = sqrt((8.0*Ru*T)/(PI_NUMBER*Ms[iSpecies]));
    sig_s = 1E-20*(5E4*5E4)/(T*T);

    tauP[iSpecies] = 1/(sig_s*Cs*N);

    /*--- Species relaxation time ---*/
    taus[iSpecies] = tauMW[iSpecies] + tauP[iSpecies];

    /*--- Calculate vib.-el. energies ---*/
    estar[iSpecies] = variable->CalcEve(config, T, iSpecies);

    /*--- Add species contribution to residual ---*/
    val_residual[nEv] += rhos * (estar[iSpecies] -
                                 eve_i[iSpecies]) / taus[iSpecies] * Volume;
  }

  



 // std::ofstream outfile;
//
 // outfile.open("prints.txt", std::ios_base::app); // append instead of overwrite
//
 // outfile << "val_residual[nEv]=" << val_residual[nEv] << endl; 


  /*---Set source term ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_source[iVar] = val_source[iVar]+val_residual[iVar]/Volume;

  if (implicit) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

      /*--- Rename ---*/
      rhos = V_i[RHOS_INDEX+iSpecies];
      //cat: settve first
      Cvvsst[iSpecies] = variable->CalcCvve(T, config, iSpecies);

      for (iVar = 0; iVar < nVar; iVar++) {
        val_Jacobian_i[nEv][iVar] += rhos/taus[iSpecies]*(Cvvsst[iSpecies]*dTdU_i[iVar] -
                                                          Cvve_i[iSpecies]*dTvedU_i[iVar])*Volume;
      }
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian_i[nEv][iSpecies] += (estar[iSpecies]-eve_i[iSpecies])/taus[iSpecies]*Volume;
  }
}

void CSource_NEMO::ComputeAxisymmetric(su2double *val_residual,
                                       su2double *val_source,
                                       su2double **val_Jacobian,
                                       CConfig *config) {

  unsigned short iDim, iSpecies, jSpecies, iVar, jVar;
  su2double rho, rhou, rhov, rhoEve, vel2, H, yinv;

  /*--- Calculate inverse of y coordinate ---*/
  if (Coord_i[1]!= 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;

  /*--- Rename for convenience ---*/
  rho    = V_i[RHO_INDEX];
  rhou   = U_i[nSpecies];
  rhov   = U_i[nSpecies+1];
  rhoEve = U_i[nSpecies+nDim+1];
  H      = V_i[H_INDEX];
  vel2   = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    vel2 += V_i[VEL_INDEX+iDim]*V_i[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Y[iSpecies] = V_i[RHOS_INDEX+iSpecies] / rho;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    val_residual[iSpecies] = yinv*rhov*Y[iSpecies]*Volume;
  val_residual[nSpecies]   = yinv*rhov*U_i[nSpecies]/rho*Volume;
  val_residual[nSpecies+1] = yinv*rhov*U_i[nSpecies+1]/rho*Volume;
  val_residual[nSpecies+2] = yinv*rhov*H*Volume;
  val_residual[nSpecies+3] = yinv*rhov*U_i[nSpecies+nDim+1]/rho*Volume;

  /*---Set source term ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_source[iVar] = val_source[iVar]+val_residual[iVar]/Volume;

  if (implicit) {

    /*--- Initialize ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian[iVar][jVar] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
        dYdr[iSpecies][jSpecies] = 0.0;

    /*--- Calculate additional quantities ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        dYdr[iSpecies][jSpecies] += -1/rho*Ys[iSpecies];
      }
      dYdr[iSpecies][iSpecies] += 1/rho;
    }

    /*--- Populate Jacobian ---*/

    // Species density
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        val_Jacobian[iSpecies][jSpecies] = dYdr[iSpecies][jSpecies]*rhov;
      }
      val_Jacobian[iSpecies][nSpecies+1] = Y[iSpecies];
    }

    // X-momentum
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian[nSpecies][iSpecies] = -rhou*rhov/(rho*rho);
    val_Jacobian[nSpecies][nSpecies] = rhov/rho;
    val_Jacobian[nSpecies][nSpecies+1] = rhou/rho;

    // Y-momentum
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian[nSpecies+1][iSpecies] = -rhov*rhov/(rho*rho);
    val_Jacobian[nSpecies+1][nSpecies+1] = 2*rhov/rho;

    // Energy
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian[nSpecies+nDim][iSpecies]      = -H*rhov/rho + dPdU_i[iSpecies]*rhov/rho;
    val_Jacobian[nSpecies+nDim][nSpecies]        = dPdU_i[nSpecies]*rhov/rho;
    val_Jacobian[nSpecies+nDim][nSpecies+1]      = H + dPdU_i[nSpecies+1]*rhov/rho;
    val_Jacobian[nSpecies+nDim][nSpecies+nDim]   = (1+dPdU_i[nSpecies+nDim])*rhov/rho;
    val_Jacobian[nSpecies+nDim][nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1]*rhov/rho;

    // Vib-el energy
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian[nSpecies+nDim+1][iSpecies] = -rhoEve*rhov/(rho*rho);
    val_Jacobian[nSpecies+nDim+1][nSpecies+1] = rhoEve/rho;
    val_Jacobian[nSpecies+nDim+1][nSpecies+nDim+1] = rhov/rho;

    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian[iVar][jVar] *= yinv*Volume;
  }
}

