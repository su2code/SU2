/*!
 * \file CTNE2EulerVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon, S.R. Copeland, W. Maier
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/variables/CTNE2EulerVariable.hpp"
#include <math.h>

CTNE2EulerVariable::CTNE2EulerVariable(unsigned long npoint,
                                       unsigned long val_ndim,
                                       unsigned long val_nvar,
                                       unsigned long val_nprimvar,
                                       unsigned long val_nprimvargrad,
                                       CConfig *config) : CVariable(npoint,
                                                                    val_ndim,
                                                                    val_nvar,
                                                                    config) {

  nDim         = val_ndim;
  nVar         = val_nvar;

  nPrimVar     = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;

  nSpecies     = config->GetnSpecies();
  ionization   = config->GetIonization();

  /*--- Array initialization ---*/
  Limiter_Primitive.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Limiter.resize(nPoint,nVar) = su2double(0.0);
  Primitive.resize(nPoint,nPrimVar) = su2double(0.0);
  Gradient_Primitive.resize(nPoint,nPrimVarGrad,nDim,0.0);
  Gradient.resize(nPoint,nVar,nDim,0.0);

  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;

}




CTNE2EulerVariable::CTNE2EulerVariable(su2double val_pressure,
                                       su2double *val_massfrac,
                                       su2double *val_mach,
                                       su2double val_temperature,
                                       su2double val_temperature_ve,
                                       unsigned long npoint,
                                       unsigned long ndim,
                                       unsigned long nvar,
                                       unsigned long nvarprim,
                                       unsigned long nvarprimgrad,
                                       CConfig *config) : CVariable(npoint,
                                                                    ndim,
                                                                    nvar,
                                                                    config   ) {

  unsigned short iEl, iDim, iSpecies, iVar, nEl, nHeavy, nMGSmooth;
  unsigned short *nElStates;
  su2double *xi, *Ms, *thetav, **thetae, **g, *hf, *Tref;
  su2double rhoE, rhoEve, Ev, Ee, Ef, T, Tve, rho, rhoCvtr, rhos;
  su2double RuSI, Ru, sqvel, num, denom, conc, soundspeed;

  /*--- Setting variable amounts ---*/
  nSpecies     = config->GetnSpecies();
  nDim         = ndim;
  nPrimVar     = nvarprim;
  nPrimVarGrad = nvarprimgrad;

  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;

  /*--- Allocate & initialize residual vectors ---*/

  Res_TruncError.resize(nPoint,nVar) = su2double(0.0);

  /*--- Only for residual smoothing (multigrid) ---*/

  for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (config->GetMG_CorrecSmooth(iMesh) > 0) {
      Residual_Sum.resize(nPoint,nVar);
      Residual_Old.resize(nPoint,nVar);
      break;
    }
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian.resize(nPoint,nVar);

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Limiter.resize(nPoint,nVar) = su2double(0.0);

  Solution_Max.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Solution_Min.resize(nPoint,nPrimVarGrad) = su2double(0.0);

  /*--- Primitive and secondary variables ---*/

  Primitive.resize(nPoint,nPrimVar) = su2double(0.0);

  dPdU.resize(nPoint, nVar)      = su2double(0.0);
  dTdU.resize(nPoint, nVar)      = su2double(0.0);
  dTvedU.resize(nPoint, nVar)    = su2double(0.0);
  Cvves.resize(nPoint, nSpecies) = su2double(0.0);
  eves.resize(nPoint, nSpecies)  = su2double(0.0);
  //Secondary.resize(nPoint,nSecondaryVar) = su2double(0.0);
  Source.resize(nPoint,nVar) = su2double(0.0);
  /*--- Compressible flow, gradients primitive variables ---*/

  Gradient_Primitive.resize(nPoint,nPrimVarGrad,nDim,0.0);
  Gradient.resize(nPoint,nVar,nDim,0.0);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
   Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  //if (config->GetMultizone_Problem())
  //  Set_BGSSolution_k();

  Velocity2.resize(nPoint) = su2double(0.0);
  Max_Lambda_Inv.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);
  Lambda.resize(nPoint) = su2double(0.0);
  Sensor.resize(nPoint) = su2double(0.0);

  /* Under-relaxation parameter. */
  //UnderRelaxation.resize(nPoint) = su2double(1.0);
  //LocalCFL.resize(nPoint) = su2double(0.0);

  /* Non-physical point (first-order) initialization. */
  Non_Physical.resize(nPoint) = false;
  //Non_Physical_Counter.resize(nPoint) = 0;

  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Load variables from the config class --*/
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
  g         = config->GetElDegeneracy();       // Degeneracy of electron states
  nElStates = config->GetnElStates();          // Number of electron states
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]

  /*--- Rename & initialize for convenience ---*/
  RuSI      = UNIVERSAL_GAS_CONSTANT;          // Universal gas constant [J/(mol*K)]
  Ru        = 1000.0*RuSI;                     // Universal gas constant [J/(kmol*K)]
  Tve       = val_temperature_ve;              // Vibrational temperature [K]
  T         = val_temperature;                 // Translational-rotational temperature [K]

  /*--- Loop over all points --*/
  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint){

    /*--- Reset values to zero ---*/
    sqvel     = 0.0;                             // Velocity^2 [m2/s2]
    rhoE      = 0.0;                             // Mixture total energy per mass [J/kg]
    rhoEve    = 0.0;                             // Mixture vib-el energy per mass [J/kg]
    rhoCvtr   = 0.0;                             // Mixture Cv (trans-rot) per mass
    denom     = 0.0;
    conc      = 0.0;

    /*--- Calculate mixture density from supplied primitive quantities ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
      denom += val_massfrac[iSpecies] * (Ru/Ms[iSpecies]) * T;

    for (iSpecies = 0; iSpecies < nEl; iSpecies++)
      denom += val_massfrac[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;

    rho = val_pressure / denom;

    /*--- Calculate sound speed and extract velocities ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      conc += val_massfrac[iSpecies]*rho/Ms[iSpecies];
      rhoCvtr += rho*val_massfrac[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    }

    soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * val_pressure/rho);

    for (iDim = 0; iDim < nDim; iDim++)
      sqvel += val_mach[iDim]*soundspeed * val_mach[iDim]*soundspeed;

    /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      // Species density
      rhos = val_massfrac[iSpecies]*rho;

      // Species formation energy
      Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

      // Species vibrational energy
      if (thetav[iSpecies] != 0.0)
        Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
      else
        Ev = 0.0;

      // Species electronic energy
      num = 0.0;
      denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);

      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
        denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      }

      Ee = Ru/Ms[iSpecies] * (num/denom);

      // Mixture total energy
      rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                      + Ev + Ee + Ef + 0.5*sqvel);

      // Mixture vibrational-electronic energy
      rhoEve += rhos * (Ev + Ee);
    }

    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
      // Species formation energy
      Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

      // Electron t-r mode contributes to mixture vib-el energy
      rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
    }

    /*--- Initialize Solution & Solution_Old vectors ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Solution(iPoint,iSpecies)     = rho*val_massfrac[iSpecies];
      Solution_Old(iPoint,iSpecies) = rho*val_massfrac[iSpecies];
    }
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution(iPoint,nSpecies+iDim)     = rho*val_mach[iDim]*soundspeed;
      Solution_Old(iPoint,nSpecies+iDim) = rho*val_mach[iDim]*soundspeed;
    }
    Solution(iPoint,nSpecies+nDim)       = rhoE;
    Solution_Old(iPoint,nSpecies+nDim)   = rhoE;
    Solution(iPoint,nSpecies+nDim+1)     = rhoEve;
    Solution_Old(iPoint,nSpecies+nDim+1) = rhoEve;

    /*--- Assign primitive variables ---*/
    Primitive(iPoint,T_INDEX)   = val_temperature;
    Primitive(iPoint,TVE_INDEX) = val_temperature_ve;
    Primitive(iPoint,P_INDEX)   = val_pressure;
  }

  if (config->GetError_Estimate() || config->GetKind_SU2() == SU2_MET) {
    AnisoMetr.resize(nPoint,3*(nDim-1)) = su2double(0.0);
    if(config->GetViscous()) {
      AnisoViscGrad.resize(nPoint,nDim*nVar*nDim) = su2double(0.0);
      AnisoViscHess.resize(nPoint,3*(nDim-1)*nVar*nDim) = su2double(0.0);
    }
    if(config->GetAdap_Source()) {
      AnisoSourceGrad.resize(nPoint,nDim*nVar) = su2double(0.0);
      AnisoSourceHess.resize(nPoint, 3*(nDim-1)*nVar) = su2double(0.0);
    }
  }

}

void CTNE2EulerVariable::SetGradient_PrimitiveZero() {
  Gradient_Primitive.storage.setConstant(0.0);
}

bool CTNE2EulerVariable::SetDensity(unsigned long iPoint) {

  unsigned short iSpecies;
  su2double Density;

  Density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Primitive(iPoint,RHOS_INDEX+iSpecies) = Solution(iPoint,iSpecies);
    Density += Solution(iPoint,iSpecies);
  }
  Primitive(iPoint,RHO_INDEX) = Density;

  return true;
}

void CTNE2EulerVariable::SetVelocity2(unsigned long iPoint) {

  unsigned short iDim;

  Velocity2(iPoint) = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Primitive(iPoint,VEL_INDEX+iDim) = Solution(iPoint,nSpecies+iDim) / Primitive(iPoint,RHO_INDEX);
    Velocity2(iPoint) +=  Solution(iPoint,nSpecies+iDim)*Solution(iPoint,nSpecies+iDim)
        / (Primitive(iPoint,RHO_INDEX)*Primitive(iPoint,RHO_INDEX));
  }
}

bool CTNE2EulerVariable::SetTemperature(unsigned long iPoint, CConfig *config) {

  // Note: Requires previous call to SetDensity()
  unsigned short iEl, iSpecies, iDim, nHeavy, nEl, iIter, maxIter, *nElStates;
  su2double *xi, *Ms, *thetav, **thetae, **g, *hf, *Tref;
  su2double rho, rhoE, rhoEve, rhoEve_t, rhoE_ref, rhoE_f;
  su2double evs, eels;
  su2double RuSI, Ru, sqvel, rhoCvtr, rhoCvve;
  su2double Cvvs, Cves, Tve, Tve2, Tve_o;
  su2double f, df, tol;
  su2double exptv, thsqr, thoTve;
  su2double num, denom, num2, num3;
  su2double Tmin, Tmax, Tvemin, Tvemax;

  /*--- Set tolerance for Newton-Raphson method ---*/
  tol     = 1.0E-4;
  maxIter = 100;

  /*--- Set temperature clipping values ---*/
  Tmin   = 100.0;
  Tmax   = 6E4;
  Tvemin = 100.0;
  Tvemax = 4E4;

  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Load variables from the config class --*/
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  nElStates = config->GetnElStates();

  /*--- Rename & initialize for convenience ---*/
  RuSI     = UNIVERSAL_GAS_CONSTANT;           // Universal gas constant [J/(mol*K)]
  Ru       = 1000.0*RuSI;                      // Universal gas constant [J/(kmol*k)]
  rho      = Primitive(iPoint,RHO_INDEX);             // Mixture density [kg/m3]
  rhoE     = Solution(iPoint,nSpecies+nDim);          // Density * energy [J/m3]
  rhoEve   = Solution(iPoint,nSpecies+nDim+1);        // Density * energy_ve [J/m3]
  rhoE_f   = 0.0;                              // Density * formation energy [J/m3]
  rhoE_ref = 0.0;                              // Density * reference energy [J/m3]
  rhoCvtr  = 0.0;                              // Mix spec. heat @ const. volume [J/(kg*K)]
  sqvel    = 0.0;                              // Velocity^2 [m2/s2]

  /*--- Error checking ---*/
  if (rhoE < 0.0)
    rhoE = EPS;
  if (rhoEve < 0.0)
    rhoEve = EPS;

  /*--- Calculate mixture properties (heavy particles only) ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr  += Solution(iPoint,iSpecies) * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    rhoE_ref += Solution(iPoint,iSpecies) * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
    rhoE_f   += Solution(iPoint,iSpecies) * (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
  }
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel    += (Solution(iPoint,nSpecies+iDim)/rho) * (Solution(iPoint,nSpecies+iDim)/rho);

  /*--- Calculate translational-rotational temperature ---*/
  Primitive(iPoint,T_INDEX) = (rhoE - rhoEve - rhoE_f + rhoE_ref - 0.5*rho*sqvel) / rhoCvtr;

  /*--- Calculate vibrational-electronic temperature ---*/
  // NOTE: Cannot write an expression explicitly in terms of Tve
  // NOTE: We use Newton-Raphson to iteratively find the value of Tve
  // NOTE: Use T as an initial guess
  Tve   = Primitive(iPoint,TVE_INDEX);
  Tve_o = Primitive(iPoint,TVE_INDEX);

  for (iIter = 0; iIter < maxIter; iIter++) {
    rhoEve_t = 0.0;
    rhoCvve  = 0.0;
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {

      /*--- Vibrational energy ---*/
      if (thetav[iSpecies] != 0.0) {

        /*--- Rename for convenience ---*/
        thoTve = thetav[iSpecies]/Tve;
        exptv = exp(thetav[iSpecies]/Tve);
        thsqr = thetav[iSpecies]*thetav[iSpecies];

        /*--- Calculate vibrational energy ---*/
        evs  = Ru/Ms[iSpecies] * thetav[iSpecies] / (exptv - 1.0);

        /*--- Calculate species vibrational specific heats ---*/
        Cvvs  = Ru/Ms[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));

        /*--- Add contribution ---*/
        rhoEve_t += Solution(iPoint,iSpecies) * evs;
        rhoCvve  += Solution(iPoint,iSpecies) * Cvvs;
      }
      /*--- Electronic energy ---*/
      if (nElStates[iSpecies] != 0) {
        num = 0.0; num2 = 0.0;
        denom = g[iSpecies][0] * exp(-thetae[iSpecies][0]/Tve);
        num3  = g[iSpecies][0] * (thetae[iSpecies][0]/(Tve*Tve))*exp(-thetae[iSpecies][0]/Tve);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          thoTve = thetae[iSpecies][iEl]/Tve;
          exptv = exp(-thetae[iSpecies][iEl]/Tve);

          num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exptv;
          denom += g[iSpecies][iEl] * exptv;
          num2  += g[iSpecies][iEl] * (thoTve*thoTve) * exptv;
          num3  += g[iSpecies][iEl] * thoTve/Tve * exptv;
        }
        eels = Ru/Ms[iSpecies] * (num/denom);
        Cves = Ru/Ms[iSpecies] * (num2/denom - num*num3/(denom*denom));

        rhoEve_t += Solution(iPoint,iSpecies) * eels;
        rhoCvve  += Solution(iPoint,iSpecies) * Cves;
      }
    }
    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
      Cves = 3.0/2.0 * Ru/Ms[nSpecies-1];
      rhoEve_t += Solution(iPoint,nSpecies-1) * Cves * Tve;
      rhoCvve += Solution(iPoint,nSpecies-1) * Cves;
    }

    /*--- Determine function f(Tve) and df/dTve ---*/
    f  = rhoEve - rhoEve_t;
    df = -rhoCvve;
    Tve2 = Tve - (f/df)*0.5;

    /*--- Check for non-physical conditions ---*/
    if ((Tve2 != Tve2) || (Tve2 < 0))
      Tve2 = 1.4*Primitive(iPoint,T_INDEX);

    /*--- Check for convergence ---*/
    if (fabs(Tve2-Tve) < tol) break;
    if (iIter == maxIter-1) {
      cout << "WARNING!!! Tve convergence not reached!" << endl;
      cout << "rhoE: " << rhoE << endl;
      cout << "rhoEve: " << rhoEve << endl;
      cout << "T: " << Primitive(iPoint,T_INDEX) << endl;
      cout << "Tve2: " << Tve2 << endl;
      cout << "Tve_o: " << Tve_o << endl;
      Tve2 = Tve_o;
      break;
    }
    Tve = Tve2;
  }

  Primitive(iPoint,TVE_INDEX) = Tve2;

  /*--- Error checking ---*/
  if (Primitive(iPoint,T_INDEX) <= Tmin) {
    cout << "WARNING: T = " << Primitive(iPoint,T_INDEX) << "\t -- Clipping at: " << Tmin << endl;
    Primitive(iPoint,T_INDEX) = Tmin;
  } else if (Primitive(iPoint,T_INDEX) >= Tmax) {
    cout << "WARNING: T = " << Primitive(iPoint,T_INDEX) << "\t -- Clipping at: " << Tmax << endl;
    Primitive(iPoint,T_INDEX) = Tmax;
  }
  if (Primitive(iPoint,TVE_INDEX) <= Tvemin) {
    cout << "WARNING: Tve = " << Primitive(iPoint,TVE_INDEX) << "\t -- Clipping at: " << Tvemin << endl;
    Primitive(iPoint,TVE_INDEX) = Tvemin;
  } else if (Primitive(iPoint,TVE_INDEX) >= Tvemax) {
    cout << "WARNING: Tve = " << Primitive(iPoint,TVE_INDEX) << "\t -- Clipping at: " << Tvemax << endl;
    Primitive(iPoint,TVE_INDEX) = Tvemax;
  }

  /*--- Assign Gas Properties ---*/
  Primitive(iPoint,RHOCVTR_INDEX) = rhoCvtr;
  Primitive(iPoint,RHOCVVE_INDEX) = rhoCvve;

  /*--- Check that the solution is physical ---*/
  if ((Primitive(iPoint,T_INDEX) > 0.0) && (Primitive(iPoint,TVE_INDEX) > 0.0)) return false;
  else return true;
}

bool CTNE2EulerVariable::SetPressure(unsigned long iPoint, CConfig *config) {

  // NOTE: Requires computation of trans-rot & vib-el temperatures.

  unsigned short iSpecies, nHeavy, nEl;
  su2double *Ms;
  su2double P, RuSI, Ru;

  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Read gas mixture properties from config ---*/
  Ms = config->GetMolar_Mass();

  /*--- Rename for convenience ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;

  /*--- Solve for mixture pressure using ideal gas law & Dalton's law ---*/
  // Note: If free electrons are present, use Tve for their partial pressure
  P = 0.0;
  for(iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    P += Solution(iPoint,iSpecies) * Ru/Ms[iSpecies] * Primitive(iPoint,T_INDEX);

  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    P += Solution(iPoint,nSpecies-1) * Ru/Ms[nSpecies-1] * Primitive(iPoint,TVE_INDEX);

  /*--- Store computed values and check for a physical solution ---*/
  Primitive(iPoint,P_INDEX) = P;
  if (Primitive(iPoint,P_INDEX) > 0.0) return false;
  else return true;
}

bool CTNE2EulerVariable::SetSoundSpeed(unsigned long iPoint) {

  unsigned short iSpecies, iDim;
  su2double radical2;

  radical2 = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    radical2 += Primitive(iPoint,RHOS_INDEX+iSpecies)/Primitive(iPoint, RHO_INDEX) * dPdU(iPoint,iSpecies);
  for (iDim = 0; iDim < nDim; iDim++)
    radical2 += Primitive(iPoint,VEL_INDEX+iDim)*dPdU(iPoint,nSpecies+iDim);
  radical2 += (Solution(iPoint,nSpecies+nDim)+Primitive(iPoint,P_INDEX))/Primitive(iPoint,RHO_INDEX) * dPdU(iPoint,nSpecies+nDim);
  radical2 += Solution(iPoint,nSpecies+nDim+1)/Primitive(iPoint,RHO_INDEX) * dPdU(iPoint,nSpecies+nDim+1);

  if (radical2 < 0.0) return true;
  else {Primitive(iPoint,A_INDEX) = sqrt(radical2); return false; }

}

void CTNE2EulerVariable::CalcdPdU(su2double *V, su2double *val_eves,
                                  CConfig *config, su2double *val_dPdU) {

  // Note: Requires SetDensity(), SetTemperature(), SetPressure(), & SetGasProperties()
  // Note: Electron energy not included properly.

  unsigned short iDim, iSpecies, iEl, nHeavy, nEl, *nElStates;
  su2double *Ms, *Tref, *hf, *xi, *thetav, **thetae, **g;
  su2double RuSI, Ru, RuBAR, CvtrBAR, rhoCvtr, rhoCvve, Cvtrs, rho_el, sqvel, conc;
  su2double rho, rhos, T, Tve, ef;
  su2double num, denom;

  if (val_dPdU == NULL) {
    cout << "ERROR: CalcdPdU - Array dPdU not allocated!" << endl;
    exit(1);
  }

  /*--- Determine the number of heavy species ---*/
  if (ionization) {
    nHeavy = nSpecies-1;
    nEl    = 1;
    rho_el = V[RHOS_INDEX+nSpecies-1];
  } else {
    nHeavy = nSpecies;
    nEl    = 0;
    rho_el = 0.0;
  }

  /*--- Read gas mixture properties from config ---*/
  Ms        = config->GetMolar_Mass();
  Tref      = config->GetRefTemperature();
  hf        = config->GetEnthalpy_Formation();
  xi        = config->GetRotationModes();
  thetav    = config->GetCharVibTemp();
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  nElStates = config->GetnElStates();

  /*--- Rename for convenience ---*/
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;
  T       = V[T_INDEX];
  Tve     = V[TVE_INDEX];
  rho     = V[RHO_INDEX];
  rhoCvtr = V[RHOCVTR_INDEX];
  rhoCvve = V[RHOCVVE_INDEX];

  /*--- Pre-compute useful quantities ---*/
  RuBAR   = 0.0;
  CvtrBAR = 0.0;
  sqvel   = 0.0;
  conc    = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim] * V[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    CvtrBAR += V[RHOS_INDEX+iSpecies]*(3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    conc    += V[RHOS_INDEX+iSpecies]/Ms[iSpecies];
  }

  // Species density
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhos  = V[RHOS_INDEX+iSpecies];
    ef    = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    Cvtrs = (3.0/2.0+xi[iSpecies]/2.0)*Ru/Ms[iSpecies];

    val_dPdU[iSpecies] =  T*Ru/Ms[iSpecies] + Ru*conc/rhoCvtr *
        (-Cvtrs*(T-Tref[iSpecies]) -
         ef + 0.5*sqvel);
  }
  if (ionization) {
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      //      evibs = Ru/Ms[iSpecies] * thetav[iSpecies]/(exp(thetav[iSpecies]/Tve)-1.0);
      //      num = 0.0;
      //      denom = g[iSpecies][0] * exp(-thetae[iSpecies][0]/Tve);
      //      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      //        num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      //        denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      //      }
      //      eels = Ru/Ms[iSpecies] * (num/denom);

      val_dPdU[iSpecies] -= rho_el * Ru/Ms[nSpecies-1] * (val_eves[iSpecies])/rhoCvve;
    }
    ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1]*Tref[nSpecies-1];
    val_dPdU[nSpecies-1] = Ru*conc/rhoCvtr * (-ef + 0.5*sqvel)
        + Ru/Ms[nSpecies-1]*Tve
        - rho_el*Ru/Ms[nSpecies-1] * (-3.0/2.0*Ru/Ms[nSpecies-1]*Tve)/rhoCvve;
  }
  // Momentum
  for (iDim = 0; iDim < nDim; iDim++)
    val_dPdU[nSpecies+iDim] = -conc*Ru*V[VEL_INDEX+iDim]/rhoCvtr;

  // Total energy
  val_dPdU[nSpecies+nDim]   = conc*Ru / rhoCvtr;

  // Vib.-el energy
  val_dPdU[nSpecies+nDim+1] = -val_dPdU[nSpecies+nDim]
      + rho_el*Ru/Ms[nSpecies-1]*1.0/rhoCvve;
}

su2double CTNE2EulerVariable::CalcEve(CConfig *config, su2double val_Tve,
                                      unsigned short val_Species) {

  unsigned short iEl, *nElStates;
  su2double *Ms, *thetav, **thetae, **g, *hf, *Tref, RuSI, Ru;
  su2double Tve, Ev, Eel, Ef;
  su2double num, denom;

  /*--- Read gas mixture properties from config ---*/
  Ms        = config->GetMolar_Mass();

  /*--- Rename for convenience ---*/
  RuSI  = UNIVERSAL_GAS_CONSTANT;
  Ru    = 1000.0*RuSI;
  Tve   = val_Tve;

  /*--- Electron species energy ---*/
  if ((ionization) && (val_Species == nSpecies-1)) {

    /*--- Get quantities from CConfig ---*/
    Tref = config->GetRefTemperature();
    hf   = config->GetEnthalpy_Formation();

    /*--- Calculate formation energy ---*/
    Ef = hf[val_Species] - Ru/Ms[val_Species] * Tref[val_Species];

    /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
    Eel = (3.0/2.0) * Ru/Ms[val_Species] * (Tve - Tref[val_Species]) + Ef;
    Ev  = 0.0;

  }

  /*--- Heavy particle energy ---*/
  else {

    /*--- Read from CConfig ---*/
    thetav    = config->GetCharVibTemp();
    thetae    = config->GetCharElTemp();
    g         = config->GetElDegeneracy();
    nElStates = config->GetnElStates();

    /*--- Calculate vibrational energy (harmonic-oscillator model) ---*/
    if (thetav[val_Species] != 0.0)
      Ev = Ru/Ms[val_Species] * thetav[val_Species] / (exp(thetav[val_Species]/Tve)-1.0);
    else
      Ev = 0.0;

    /*--- Calculate electronic energy ---*/
    num = 0.0;
    denom = g[val_Species][0] * exp(-thetae[val_Species][0]/Tve);
    for (iEl = 1; iEl < nElStates[val_Species]; iEl++) {
      num   += g[val_Species][iEl] * thetae[val_Species][iEl] * exp(-thetae[val_Species][iEl]/Tve);
      denom += g[val_Species][iEl] * exp(-thetae[val_Species][iEl]/Tve);
    }
    Eel = Ru/Ms[val_Species] * (num/denom);
  }

  return Ev + Eel;
}

su2double CTNE2EulerVariable::CalcHs(CConfig *config, su2double val_T,
                                     su2double val_eves, unsigned short val_Species) {

  su2double RuSI, Ru, *xi, *Ms, *hf, *Tref, T, eve, ef, hs;

  /*--- Read from config ---*/
  xi   = config->GetRotationModes();
  Ms   = config->GetMolar_Mass();
  hf   = config->GetEnthalpy_Formation();
  Tref = config->GetRefTemperature();

  /*--- Rename for convenience ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  T    = val_T;

  /*--- Calculate vibrational-electronic energy per unit mass ---*/
  eve = val_eves;

  /*--- Calculate formation energy ---*/
  ef = hf[val_Species] - Ru/Ms[val_Species]*Tref[val_Species];

  hs = Ru/Ms[val_Species]*T
      + (3.0/2.0+xi[val_Species]/2.0)*Ru/Ms[val_Species]*T
      + hf[val_Species] + eve;

  return hs;
}

su2double CTNE2EulerVariable::CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species) {

  unsigned short iEl, *nElStates;
  su2double *Ms, *thetav, **thetae, **g, RuSI, Ru;
  su2double thoTve, exptv, thsqr, Cvvs, Cves;
  su2double Tve;
  su2double num, num2, num3, denom;

  /*--- Read from config ---*/
  thetav    = config->GetCharVibTemp();
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  Ms        = config->GetMolar_Mass();
  nElStates = config->GetnElStates();

  /*--- Rename for convenience ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Tve  = val_Tve;

  /*--- If requesting electron specific heat ---*/
  if (ionization && val_Species == nSpecies-1) {
    Cvvs = 0.0;
    Cves = 3.0/2.0 * Ru/Ms[nSpecies-1];
  }

  /*--- Heavy particle specific heat ---*/
  else {

    /*--- Vibrational energy ---*/
    if (thetav[val_Species] != 0.0) {
      thoTve = thetav[val_Species]/Tve;
      exptv = exp(thetav[val_Species]/Tve);
      thsqr = thetav[val_Species]*thetav[val_Species];
      Cvvs  = Ru/Ms[val_Species] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
    } else {
      Cvvs = 0.0;
    }

    /*--- Electronic energy ---*/
    if (nElStates[val_Species] != 0) {
      num = 0.0; num2 = 0.0;
      denom = g[val_Species][0] * exp(-thetae[val_Species][0]/Tve);
      num3  = g[val_Species][0] * (thetae[val_Species][0]/(Tve*Tve))*exp(-thetae[val_Species][0]/Tve);
      for (iEl = 1; iEl < nElStates[val_Species]; iEl++) {
        thoTve = thetae[val_Species][iEl]/Tve;
        exptv = exp(-thetae[val_Species][iEl]/Tve);

        num   += g[val_Species][iEl] * thetae[val_Species][iEl] * exptv;
        denom += g[val_Species][iEl] * exptv;
        num2  += g[val_Species][iEl] * (thoTve*thoTve) * exptv;
        num3  += g[val_Species][iEl] * thoTve/Tve * exptv;
      }
      Cves = Ru/Ms[val_Species] * (num2/denom - num*num3/(denom*denom));
    } else {
      Cves = 0.0;
    }

  }
  return Cvvs + Cves;
}

void CTNE2EulerVariable::CalcdTdU(su2double *V, CConfig *config,
                                  su2double *val_dTdU) {

  unsigned short iDim, iSpecies, nHeavy, nEl;
  su2double *Ms, *xi, *Tref, *hf;
  su2double v2, ef, T, Cvtrs, rhoCvtr, RuSI, Ru;

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Get gas properties from config settings ---*/
  Ms   = config->GetMolar_Mass();
  xi   = config->GetRotationModes();
  hf   = config->GetEnthalpy_Formation();
  Tref = config->GetRefTemperature();

  /*--- Rename for convenience ---*/
  rhoCvtr = V[RHOCVTR_INDEX];
  T       = V[T_INDEX];
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;

  /*--- Calculate supporting quantities ---*/
  v2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    v2 += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];

  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    ef    = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    Cvtrs = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    val_dTdU[iSpecies]   = (-ef + 0.5*v2 + Cvtrs*(Tref[iSpecies]-T)) / rhoCvtr;
  }
  if (ionization) {
    cout << "CTNE2Variable: NEED TO IMPLEMENT dTdU for IONIZED MIX" << endl;
    exit(1);
  }

  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTdU[nSpecies+iDim] = -V[VEL_INDEX+iDim] / V[RHOCVTR_INDEX];

  /*--- Energy derivatives ---*/
  val_dTdU[nSpecies+nDim]   =  1.0 / V[RHOCVTR_INDEX];
  val_dTdU[nSpecies+nDim+1] = -1.0 / V[RHOCVTR_INDEX];

}

void CTNE2EulerVariable::CalcdTvedU(su2double *V, su2double *val_eves, CConfig *config,
                                    su2double *val_dTvedU) {

  unsigned short iDim, iSpecies;
  su2double rhoCvve;

  /*--- Rename for convenience ---*/
  rhoCvve = V[RHOCVVE_INDEX];

  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    val_dTvedU[iSpecies] = -val_eves[iSpecies]/rhoCvve;
  }
  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTvedU[nSpecies+iDim] = 0.0;

  /*--- Energy derivatives ---*/
  val_dTvedU[nSpecies+nDim]   = 0.0;
  val_dTvedU[nSpecies+nDim+1] = 1.0 / rhoCvve;

}

bool CTNE2EulerVariable::SetPrimVar_Compressible(unsigned long iPoint, CConfig *config) {

  bool nonPhys, bkup;
  unsigned short iVar;

  /*--- Convert conserved to primitive variables ---*/
  nonPhys = Cons2PrimVar(config, Solution[iPoint], Primitive[iPoint],
                         dPdU[iPoint], dTdU[iPoint], dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);
  if (nonPhys) {
    // cout << "Non-physical at " << iPoint << "." << endl;
    for (iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);
    bkup = Cons2PrimVar(config, Solution[iPoint], Primitive[iPoint], dPdU[iPoint], dTdU[iPoint],
                        dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);
    if(bkup) Prim2ConsVar(config, iPoint, Primitive[iPoint], Solution[iPoint]);
  }

  SetVelocity2(iPoint);

  return nonPhys;

}

bool CTNE2EulerVariable::Cons2PrimVar(CConfig *config, su2double *U, su2double *V,
                                      su2double *val_dPdU, su2double *val_dTdU,
                                      su2double *val_dTvedU, su2double *val_eves,
                                      su2double *val_Cvves) {

  bool ionization, errT, errTve, NRconvg, Bconvg, nonPhys;
  unsigned short iDim, iEl, iSpecies, nHeavy, nEl, iIter, maxBIter, maxNIter;
  su2double rho, rhoE, rhoEve, rhoE_f, rhoE_ref, rhoEve_min, rhoEve_max, rhoEve_t;
  su2double RuSI, Ru, sqvel, rhoCvtr, rhoCvve;
  su2double Tve, Tve2, Tve_o;
  su2double f, df, NRtol, Btol, scale;
  su2double Tmin, Tmax, Tvemin, Tvemax;
  su2double radical2;
  su2double *xi, *Ms, *hf, *Tref;

  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  // V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T

  /*--- Set booleans ---*/
  errT    = false;
  errTve  = false;
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0; Tmax   = 8E4;
  Tvemin = 50.0; Tvemax = 8E4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0E-4;      // Tolerance for the Newton-Raphson method
  Btol     = 1.0E-2;      // Tolerance for the Bisection method
  maxNIter = 999;         // Maximum Newton-Raphson iterations
  maxBIter = 999;         // Maximum Bisection method iterations
  scale    = 0.5;         // Scaling factor for Newton-Raphson step

  /*--- Read parameters from config ---*/
  xi         = config->GetRotationModes();      // Rotational modes of energy storage
  Ms         = config->GetMolar_Mass();         // Species molar mass
  Tref       = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf         = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  ionization = config->GetIonization();         // Molecule Ionization

  /*--- Rename variables for convenience ---*/
  RuSI   = UNIVERSAL_GAS_CONSTANT;    // Universal gas constant [J/(mol*K)]
  Ru     = 1000.0*RuSI;               // Universal gas constant [J/(kmol*K)]
  rhoE   = U[nSpecies+nDim];          // Density * energy [J/m3]
  rhoEve = U[nSpecies+nDim+1];        // Density * energy_ve [J/m3]

  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Assign species & mixture density ---*/
  // Note: if any species densities are < 0, these values are re-assigned
  //       in the primitive AND conserved vectors to ensure positive density
  V[RHO_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (U[iSpecies] < 0.0) {
      V[RHOS_INDEX+iSpecies] = 1E-20;
      U[iSpecies]            = 1E-20;
      nonPhys                = true;
    } else
      V[RHOS_INDEX+iSpecies] = U[iSpecies];
    V[RHO_INDEX]            += U[iSpecies];
  }

  /*--- Assign mixture velocity ---*/
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    V[VEL_INDEX+iDim] = U[nSpecies+iDim]/V[RHO_INDEX];
    sqvel            += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  }

  /*--- Translational-Rotational Temperature ---*/

  // Rename for convenience
  rho = V[RHO_INDEX];

  // Determine properties of the mixture at the current state
  rhoE_f   = 0.0;
  rhoE_ref = 0.0;
  rhoCvtr  = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr  += U[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    rhoE_ref += U[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
    rhoE_f   += U[iSpecies] * (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
  }

  // Calculate T-R temperature
  V[T_INDEX] = (rhoE - rhoEve - rhoE_f + rhoE_ref - 0.5*rho*sqvel) / rhoCvtr;
  V[RHOCVTR_INDEX] = rhoCvtr;

  // Determine if the temperature lies within the acceptable range
  if (V[T_INDEX] < Tmin) {
    V[T_INDEX] = Tmin;
    nonPhys = true;
    errT    = true;
  } else if (V[T_INDEX] > Tmax){
    V[T_INDEX] = Tmax;
    nonPhys = true;
    errT    = true;
  }


  /*--- Vibrational-Electronic Temperature ---*/

  // Check for non-physical solutions
  rhoEve_min = 0.0;
  rhoEve_max = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhoEve_min += U[iSpecies]*CalcEve(config, Tvemin, iSpecies);
    rhoEve_max += U[iSpecies]*CalcEve(config, Tvemax, iSpecies);
  }
  if (rhoEve < rhoEve_min) {
    errTve       = true;
    nonPhys      = true;
    V[TVE_INDEX] = Tvemin;
  } else if (rhoEve > rhoEve_max) {
    errTve       = true;
    nonPhys      = true;
    V[TVE_INDEX] = Tvemax;
  } else {

    AD_BEGIN_PASSIVE

    /*--- Execute a Newton-Raphson root-finding method to find Tve ---*/
    // Initialize to the translational-rotational temperature
    Tve   = V[T_INDEX];
    // // Initialize to previous Tve
    // Tve   = V[TVE_INDEX];

    // Execute the root-finding method
    NRconvg = false;
       for (iIter = 0; iIter < maxNIter; iIter++) {
         rhoEve_t = 0.0;
         rhoCvve  = 0.0;
         for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
           val_eves[iSpecies]  = CalcEve(config, Tve, iSpecies);
           val_Cvves[iSpecies] = CalcCvve(Tve, config, iSpecies);
           rhoEve_t += U[iSpecies]*val_eves[iSpecies];
           rhoCvve  += U[iSpecies]*val_Cvves[iSpecies];
         }
    
         // Find the root
         f  = U[nSpecies+nDim+1] - rhoEve_t;
         df = -rhoCvve;
         Tve2 = Tve - (f/df)*scale;
    
         // Check for nonphysical steps
         if ((Tve2 < Tvemin) || (Tve2 > Tvemax))
           break;
    //      if (Tve2 < Tvemin)
    //        Tve2 = Tvemin;
    //      else if (Tve2 > Tvemax)
    //        Tve2 = Tvemax;
    
         // Check for convergence
         if (fabs(f)/U[nSpecies+nDim+1] < NRtol) {
           NRconvg = true;
           break;
         } else {
           Tve = Tve2;
         }
       }

    // If the Newton-Raphson method has converged, assign the value of Tve.
    // Otherwise, execute a bisection root-finding method
    if (NRconvg){

      AD_END_PASSIVE

      // Recompute Eve and Cvve after search
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        val_eves[iSpecies]  = CalcEve(config, Tve, iSpecies);
        val_Cvves[iSpecies] = CalcCvve(Tve, config, iSpecies);
      }
      V[TVE_INDEX] = Tve;
    } else {

      // Assign the bounds
      Tve_o = max(Tvemin, V[T_INDEX]); // Tve should be geq T
      Tve2  = Tvemax;
      // if (rhoEve_t > rhoEve) {
      //   Tve2  = Tve;
      //   Tve_o = Tvemin;
      // }
      // else{
      //   Tve_o = Tve;
      //   Tve2  = Tvemax;
      // }

      // Execute the root-finding method
      Bconvg = false;
      for (iIter = 0; iIter < maxBIter; iIter++) {
        Tve      = (Tve_o+Tve2)/2.0;
        rhoEve_t = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies] = CalcEve(config, Tve, iSpecies);
          rhoEve_t          += U[iSpecies] * val_eves[iSpecies];
        }

        if (fabs(rhoEve_t - U[nSpecies+nDim+1])/U[nSpecies+nDim+1] < Btol) {
          V[TVE_INDEX] = Tve;
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
            val_Cvves[iSpecies] = CalcCvve(Tve, config, iSpecies);
          Bconvg = true;
          break;
        } else {
          if (rhoEve_t > rhoEve) Tve2 = Tve;
          else                  Tve_o = Tve;
        }
      }

      AD_END_PASSIVE

      // Recompute Eve and Cvve after search
      if (Bconvg) {
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies] = CalcEve(config, Tve, iSpecies);
          val_Cvves[iSpecies] = CalcCvve(Tve, config, iSpecies);
        }
      }
      // If absolutely no convergence, then assign to the TR temperature
      else {
        V[TVE_INDEX] = V[T_INDEX];
        rhoCvve  = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies]  = CalcEve(config, V[TVE_INDEX], iSpecies);
          val_Cvves[iSpecies] = CalcCvve(V[TVE_INDEX], config, iSpecies);
          rhoCvve  += U[iSpecies]*val_Cvves[iSpecies];
        }
        U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
                            - rhoE_ref + 0.5*rho*sqvel;
        U[nSpecies+nDim+1] = rhoCvve*V[TVE_INDEX];
      }
    }
  }

  /*--- Set mixture rhoCvve ---*/
  rhoCvve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoCvve += U[iSpecies]*val_Cvves[iSpecies];
  V[RHOCVVE_INDEX] = rhoCvve;

  /*--- If there are clipped temperatures, correct the energy terms ---*/
  //  if (errT) {
  //    U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
  //                       - rhoE_ref + 0.5*rho*sqvel;
  //  }
  //  if (errTve) {
  //    U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
  //                       - rhoE_ref + 0.5*rho*sqvel;
  //    U[nSpecies+nDim+1] = rhoCvve*V[TVE_INDEX];
  //  }

  /*--- Pressure ---*/
  V[P_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    V[P_INDEX] += U[iSpecies] * Ru/Ms[iSpecies] * V[T_INDEX];
  for (iEl = 0; iEl < nEl; iEl++)
    V[P_INDEX] += U[nSpecies-1] * Ru/Ms[nSpecies-1] * V[TVE_INDEX];

  if (V[P_INDEX] < 0.0) {
    V[P_INDEX] = 1E-20;
    nonPhys = true;
  }

  /*--- Partial derivatives of pressure and temperature ---*/
  CalcdPdU(  V, val_eves, config, val_dPdU  );
  CalcdTdU(  V, config, val_dTdU  );
  CalcdTvedU(V, val_eves, config, val_dTvedU);


  /*--- Sound speed ---*/
  radical2 = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    radical2 += V[RHOS_INDEX+iSpecies]/V[RHO_INDEX] * val_dPdU[iSpecies];
  for (iDim = 0; iDim < nDim; iDim++)
    radical2 += V[VEL_INDEX+iDim]*val_dPdU[nSpecies+iDim];
  radical2 += (U[nSpecies+nDim]+V[P_INDEX])/V[RHO_INDEX] * val_dPdU[nSpecies+nDim];
  radical2 += U[nSpecies+nDim+1]/V[RHO_INDEX] * val_dPdU[nSpecies+nDim+1];
  V[A_INDEX] = sqrt(radical2);

  if (radical2 < 0.0) {
    nonPhys = true;
    V[A_INDEX] = EPS;
  }


  /*--- Enthalpy ---*/
  V[H_INDEX] = (U[nSpecies+nDim] + V[P_INDEX])/V[RHO_INDEX];

  return nonPhys;
}

void CTNE2EulerVariable::Prim2ConsVar(CConfig *config, unsigned long iPoint, su2double *V, su2double *U) {
  unsigned short iDim, iEl, iSpecies, nEl, nHeavy;
  unsigned short *nElStates;
  su2double Ru, RuSI, Tve, T, sqvel, rhoE, rhoEve, Ef, Ev, Ee, rhos, rhoCvtr, num, denom;
  su2double *thetav, *Ms, *xi, *hf, *Tref;
  su2double **thetae, **g;

  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Load variables from the config class --*/
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
  g         = config->GetElDegeneracy();       // Degeneracy of electron states
  nElStates = config->GetnElStates();          // Number of electron states
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]

  /*--- Rename & initialize for convenience ---*/
  RuSI    = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(mol*K)] (SI units)
  Ru      = 1000.0*RuSI;                    // Universal gas constant [J/(kmol*K)]
  Tve     = V[TVE_INDEX];                   // Vibrational temperature [K]
  T       = V[T_INDEX];                     // Translational-rotational temperature [K]
  sqvel   = 0.0;                            // Velocity^2 [m2/s2]
  rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
  rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
  denom   = 0.0;
  rhoCvtr = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];

  /*--- Set species density ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    U[iSpecies] = V[RHOS_INDEX+iSpecies];

  /*--- Set momentum ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    U[nSpecies+iDim] = V[RHO_INDEX]*V[VEL_INDEX+iDim];

  /*--- Set the total energy ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhos = U[iSpecies];

    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
    }
    Ee = Ru/Ms[iSpecies] * (num/denom);

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);

    // Mixture vibrational-electronic energy
    rhoEve += rhos * (Ev + Ee);
  }
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
    // Species formation energy
    Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
  }

  /*--- Set energies ---*/
  U[nSpecies+nDim]   = rhoE;
  U[nSpecies+nDim+1] = rhoEve;

  return;
}

bool CTNE2EulerVariable::GradCons2GradPrimVar(CConfig *config, su2double *U,
                                              su2double *V, su2double **GradU,
                                              su2double **GradV) {

  unsigned short iSpecies, iEl, iDim, jDim, iVar, nHeavy, *nElStates;
  su2double rho, rhoCvtr, rhoCvve, T, Tve, eve, ef, eref, RuSI, Ru;
  su2double Cvvs, Cves, dCvvs, dCves;
  su2double thoTve, exptv;
  su2double An1, Bd1, Bn1, Bn2, Bn3, Bn4, A, B;
  su2double *rhou, *Grhou2;
  su2double *xi, *Ms, *Tref, *hf, *thetav;
  su2double **thetae, **g;

  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  // V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T

  /*--- Allocate arrays ---*/
  Grhou2 = new su2double[nDim];
  rhou   = new su2double[nDim];

  /*--- Determine number of heavy-particle species ---*/
  if (config->GetIonization()) nHeavy = nSpecies-1;
  else                         nHeavy = nSpecies;

  /*--- Rename for convenience ---*/
  rho     = V[RHO_INDEX];
  rhoCvtr = V[RHOCVTR_INDEX];
  rhoCvve = V[RHOCVVE_INDEX];
  T       = V[T_INDEX];
  Tve     = V[TVE_INDEX];
  xi      = config->GetRotationModes();
  Ms      = config->GetMolar_Mass();
  Tref    = config->GetRefTemperature();
  hf      = config->GetEnthalpy_Formation();
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;
  thetav  = config->GetCharVibTemp();
  g       = config->GetElDegeneracy();
  thetae  = config->GetCharElTemp();
  nElStates = config->GetnElStates();

  for (iDim = 0; iDim < nDim; iDim++)
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      GradV[iVar][iDim] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {

    Grhou2[iDim] = 0.0;

    /*--- Species density ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      GradV[RHOS_INDEX+iSpecies][iDim] = GradU[iSpecies][iDim];
      GradV[RHO_INDEX][iDim]          += GradU[iSpecies][iDim];
    }

    /*--- Velocity ---*/
    for (jDim = 0; jDim < nDim; jDim++) {
      GradV[VEL_INDEX+jDim][iDim] = (rho*GradU[nSpecies+jDim][iDim] -
          rhou[jDim]*GradV[RHO_INDEX][iDim]) / (rho*rho);
    }

    /*--- Specific Heat (T-R) ---*/
    GradV[RHOCVTR_INDEX][iDim] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      GradV[RHOCVTR_INDEX][iDim] += GradU[iSpecies][iDim]*(3.0+xi[iSpecies])/2.0 *Ru/Ms[iSpecies];

    /*--- Temperature ---*/
    // Calculate the gradient of rho*u^2
    for (jDim = 0; jDim < nDim; jDim++) {
      Grhou2[iDim] += 2.0/rho*(rhou[jDim]*GradU[nSpecies+jDim][iDim]) -
          GradV[RHO_INDEX][iDim]/(rho*rho)*(GradU[nSpecies+jDim][iDim]*
          GradU[nSpecies+jDim][iDim]);
    }
    // Calculate baseline GradT
    GradV[T_INDEX][iDim] = 1.0/rhoCvtr*(GradU[nSpecies+nDim][iDim] -
        GradU[nSpecies+nDim+1][iDim] -
        0.5*Grhou2[iDim] -
        T*GradV[RHOCVTR_INDEX][iDim]);
    // Subtract formation/reference energies
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      eref = (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
      ef   = (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
      GradV[T_INDEX][iDim] -= 1.0/rhoCvtr*(GradU[iSpecies][iDim]*(ef+eref));
    }

    /*--- Vibrational-electronic temperature ---*/
    GradV[TVE_INDEX][iDim] = GradU[nSpecies+nDim+1][iDim]/rhoCvve;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      eve = CalcEve(config, V[TVE_INDEX], iSpecies);
      GradV[TVE_INDEX][iDim] -= GradU[iSpecies][iDim]*eve;
    }

    /*--- Pressure ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      GradV[P_INDEX][iDim] += GradU[iSpecies][iDim]*Ru/Ms[iSpecies]*T +
          U[iSpecies]*Ru/Ms[iSpecies]*GradV[T_INDEX][iDim];
    }

    /*--- Enthalpy ---*/
    GradV[H_INDEX][iDim] = rho*(GradU[nSpecies+nDim][iDim] + GradV[P_INDEX][iDim]) -
        (U[nSpecies+nDim]+V[P_INDEX])*GradV[RHO_INDEX][iDim] / (rho*rho);

    /*--- Specific Heat (V-E) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      // Vibrational energy specific heat
      if (thetav[iSpecies] != 0) {
        Cvvs = Ru/Ms[iSpecies]*(thetav[iSpecies]*thetav[iSpecies]/(Tve*Tve))*
            exp(thetav[iSpecies]/Tve)/((exp(thetav[iSpecies]/Tve)-1)*
                                       (exp(thetav[iSpecies]/Tve)-1));
        dCvvs = (-2/Tve - thetav[iSpecies]/(Tve*Tve) +
                 2*thetav[iSpecies]*exp(thetav[iSpecies]/Tve)/
                 (Tve*Tve*(exp(thetav[iSpecies]/Tve)-1)))*Cvvs;
      } else {
        Cvvs = 0.0;
        dCvvs = 0.0;
      }


      // Electronic energy specific heat
      An1 = 0.0;
      Bn1 = 0.0;
      Bn2 = g[iSpecies][0]*thetae[iSpecies][0]/(Tve*Tve)*exp(-thetae[iSpecies][0]/Tve);
      Bn3 = g[iSpecies][0]*(thetae[iSpecies][0]*thetae[iSpecies][0]/
          (Tve*Tve*Tve*Tve))*exp(-thetae[iSpecies][0]/Tve);
      Bn4 = 0.0;
      Bd1 = g[iSpecies][0]*exp(-thetae[iSpecies][0]/Tve);
      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        thoTve = thetae[iSpecies][iEl]/Tve;
        exptv = exp(-thetae[iSpecies][iEl]/Tve);

        An1 += g[iSpecies][iEl]*thoTve*thoTve*exptv;
        Bn1 += g[iSpecies][iEl]*thetae[iSpecies][iEl]*exptv;
        Bn2 += g[iSpecies][iEl]*thoTve/Tve*exptv;
        Bn3 += g[iSpecies][iEl]*thoTve*thoTve/(Tve*Tve)*exptv;
        Bn4 += g[iSpecies][iEl]*thoTve*thoTve*thoTve/Tve*exptv;
        Bd1 += g[iSpecies][iEl]*exptv;
      }
      A = An1/Bd1;
      B = Bn1*Bn2/(Bd1*Bd1);
      Cves = Ru/Ms[iSpecies]*(A-B);

      dCves = Ru/Ms[iSpecies]*(-2.0/Tve*(A-B) - 2*Bn2/Bd1*(A-B) -
                               Bn1*Bn3/(Bd1*Bd1) + Bn4/Bd1);

      GradV[RHOCVVE_INDEX][iDim] += V[RHOS_INDEX+iSpecies]*(dCvvs+dCves)*GradV[TVE_INDEX][iDim] +
          GradV[RHOS_INDEX+iSpecies][iDim]*(Cvvs+Cves);

    }
  }

  delete [] Grhou2;
  delete [] rhou;
  return false;
}

void CTNE2EulerVariable::SetPrimVar_Gradient(CConfig *config, unsigned long iPoint) {

  /*--- Use built-in method on TNE2 variable global types ---*/
  GradCons2GradPrimVar(config, Solution[iPoint], Primitive[iPoint],
                       Gradient[iPoint], Gradient_Primitive[iPoint]);

}

bool CTNE2EulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  cout << "NNNNNNNNNNNNNNOOOOOOOOOOOOOOOOOOOO DELETE ME" << endl;
}

void CTNE2EulerVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {

   /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

   SetdPdrho_e(iPoint, FluidModel->GetdPdrho_e());
   SetdPde_rho(iPoint, FluidModel->GetdPde_rho());

}

void CTNE2EulerVariable::SetSolution_New() { Solution_New = Solution; }
