/*!
 * \file CTNE2NSVariable.cpp
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

#include "../../include/variables/CTNE2NSVariable.hpp"
#include <math.h>

CTNE2NSVariable::CTNE2NSVariable(unsigned long val_ndim,
                                 unsigned long val_nvar,
                                 unsigned long val_nprimvar,
                                 unsigned long val_nprimvargrad,
                                 unsigned long  npoint,
                                 CConfig *config) : CTNE2EulerVariable(npoint,
                                                                       val_ndim,
                                                                       val_nvar,
                                                                       val_nprimvar,
                                                                       val_nprimvargrad,
                                                                       config) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();

  DiffusionCoeff.resize(nPoint,nSpecies)  = su2double(0.0);
  Dij.resize(nPoint, nSpecies, nSpecies, 0.0);
}

CTNE2NSVariable::CTNE2NSVariable(su2double val_pressure, su2double *val_massfrac,
                                 su2double *val_mach, su2double val_temperature,
                                 su2double val_temperature_ve,
                                 unsigned long npoint,
                                 unsigned long val_ndim,
                                 unsigned long val_nvar,
                                 unsigned long val_nvarprim,
                                 unsigned long val_nvarprimgrad,
                                 CConfig *config) : CTNE2EulerVariable(val_pressure,
                                                                       val_massfrac,
                                                                       val_mach,
                                                                       val_temperature,
                                                                       val_temperature_ve,
                                                                       npoint,
                                                                       val_ndim,
                                                                       val_nvar,
                                                                       val_nvarprim,
                                                                       val_nvarprimgrad,
                                                                       config) {


  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff.resize(nPoint,nSpecies)  = su2double(0.0);
  Dij.resize(nPoint, nSpecies, nSpecies, 0.0);
}

//CTNE2NSVariable::CTNE2NSVariable(su2double *val_solution, unsigned short val_ndim,
//                                 unsigned long val_nvar,
//                                 unsigned long val_nprimvar,
//                                 unsigned long val_nprimvargrad,
//                                 unsigned long npoint,
//                                 CConfig *config) : CTNE2EulerVariable(val_solution,
//                                                                       val_ndim,
//                                                                       val_nvar,
//                                                                       val_nprimvar,
//                                                                       val_nprimvargrad,
//                                                                       config) {

//  Temperature_Ref = config->GetTemperature_Ref();
//  Viscosity_Ref   = config->GetViscosity_Ref();
//  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
//  Prandtl_Lam     = config->GetPrandtl_Lam();

//  DiffusionCoeff.resize(nPoint,nSpecies)  = su2double(0.0);
//  Dij.resize(nPoint, nSpecies, nSpecies, 0.0);
//}

void CTNE2NSVariable::SetDiffusionCoeff_GuptaYos(CConfig *config) {

  unsigned short iSpecies, jSpecies, nHeavy, nEl;
  su2double rho, T, Tve, P;
  su2double *Ms, Mi, Mj, pi, RuSI, Ru, kb, gam_i, gam_j, gam_t, Theta_v;
  su2double denom, d1_ij, D_ij;
  su2double ***Omega00, Omega_ij;

  for (unsigned long iPoint =0; iPoint<nPoint; iPoint++){
    /*--- Acquire gas parameters from CConfig ---*/
    Omega00 = config->GetCollisionIntegral00();
    Ms      = config->GetMolar_Mass();
    if (ionization) {nHeavy = nSpecies-1;  nEl = 1;}
    else            {nHeavy = nSpecies;    nEl = 0;}

    /*--- Rename for convenience ---*/
    rho  = Primitive(iPoint,RHO_INDEX);
    T    = Primitive(iPoint,T_INDEX);
    Tve  = Primitive(iPoint,TVE_INDEX);
    P    = Primitive(iPoint,P_INDEX);
    pi   = PI_NUMBER;
    RuSI = UNIVERSAL_GAS_CONSTANT;
    Ru   = 1000.0*RuSI;
    kb   = BOLTZMANN_CONSTANT;

    /*--- Calculate mixture gas constant ---*/
    gam_t = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      gam_t += Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Ms[iSpecies]);
    }

    /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {

      /*--- Initialize the species diffusion coefficient ---*/
      DiffusionCoeff(iPoint,iSpecies) = 0.0;

      /*--- Calculate molar concentration ---*/
      Mi      = Ms[iSpecies];
      gam_i   = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mi);
      Theta_v = config->GetCharVibTemp(iSpecies);

      denom = 0.0;
      for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
        if (jSpecies != iSpecies) {
          Mj    = Ms[jSpecies];
          gam_j = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mj);

          /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
          Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
              * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
              + Omega00[iSpecies][jSpecies][1]*log(T)
              + Omega00[iSpecies][jSpecies][2]);

          /*--- Calculate "delta1_ij" ---*/
          d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

          /*--- Calculate heavy-particle binary diffusion coefficient ---*/
          D_ij = kb*T/(P*d1_ij);
          denom += gam_j/D_ij;
        }
      }

      if (ionization) {
        jSpecies = nSpecies-1;
        Mj       = config->GetMolar_Mass(jSpecies);
        gam_j    = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mj);

        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
            * pow(Tve, Omega00[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega00[iSpecies][jSpecies][1]*log(Tve)
            + Omega00[iSpecies][jSpecies][2]);

        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
      }

      /*--- Assign species diffusion coefficient ---*/
      DiffusionCoeff(iPoint,iSpecies) = gam_t*gam_t*Mi*(1-Mi*gam_i) / denom;
    }
    if (ionization) {
      iSpecies = nSpecies-1;
      /*--- Initialize the species diffusion coefficient ---*/
      DiffusionCoeff(iPoint,iSpecies) = 0.0;

      /*--- Calculate molar concentration ---*/
      Mi      = Ms[iSpecies];
      gam_i   = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mi);

      denom = 0.0;
      for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
        if (iSpecies != jSpecies) {
          Mj    = config->GetMolar_Mass(jSpecies);
          gam_j = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mj);

          /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
          Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
              * pow(Tve, Omega00[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
              + Omega00[iSpecies][jSpecies][1]*log(Tve)
              + Omega00[iSpecies][jSpecies][2]);

          /*--- Calculate "delta1_ij" ---*/
          d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

          /*--- Calculate heavy-particle binary diffusion coefficient ---*/
          D_ij = kb*Tve/(P*d1_ij);
          denom += gam_j/D_ij;
        }
      }
      DiffusionCoeff(iPoint,iSpecies) = gam_t*gam_t*Ms[iSpecies]*(1-Ms[iSpecies]*gam_i)
          / denom;
    }
    //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    //    DiffusionCoeff[iSpecies] = 0.0;
  }
}

void CTNE2NSVariable::SetLaminarViscosity_GuptaYos(CConfig *config) {

  unsigned short iSpecies, jSpecies, nHeavy, nEl;
  su2double rho, T, Tve;
  su2double *Ms, Mi, Mj, pi, Ru, RuSI, Na, gam_i, gam_j, denom;
  su2double ***Omega11, Omega_ij, d2_ij;

  /*--- Acquire gas parameters from CConfig ---*/
  Omega11 = config->GetCollisionIntegral11();
  Ms      = config->GetMolar_Mass();
  if (ionization) {nHeavy = nSpecies-1;  nEl = 1;}
  else            {nHeavy = nSpecies;    nEl = 0;}

  /*--- Loop over all points  ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    /*--- Rename for convenience ---*/
    rho  = Primitive(iPoint,RHO_INDEX);
    T    = Primitive(iPoint,T_INDEX);
    Tve  = Primitive(iPoint, TVE_INDEX);
    pi   = PI_NUMBER;
    RuSI = UNIVERSAL_GAS_CONSTANT;
    Ru   = 1000.0*RuSI;
    Na   = AVOGAD_CONSTANT;

    LaminarViscosity(iPoint) = 0.0;

    /*--- Mixture viscosity via Gupta-Yos approximation ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      denom = 0.0;

      /*--- Calculate molar concentration ---*/
      Mi    = Ms[iSpecies];
      gam_i = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mi);
      for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
        Mj    = Ms[jSpecies];
        gam_j = Primitive(iPoint,RHOS_INDEX+jSpecies) / (rho*Mj);

        /*--- Calculate "delta" quantities ---*/
        Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
            * pow(T, Omega11[iSpecies][jSpecies][0]*log(T)*log(T)
            + Omega11[iSpecies][jSpecies][1]*log(T)
            + Omega11[iSpecies][jSpecies][2]);
        d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

        /*--- Add to denominator of viscosity ---*/
        denom += gam_j*d2_ij;
      }
      if (ionization) {
        jSpecies = nSpecies-1;
        Mj    = Ms[jSpecies];
        gam_j = Primitive(iPoint,RHOS_INDEX+jSpecies) / (rho*Mj);

        /*--- Calculate "delta" quantities ---*/
        Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
            * pow(Tve, Omega11[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega11[iSpecies][jSpecies][1]*log(Tve)
            + Omega11[iSpecies][jSpecies][2]);
        d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

        denom += gam_j*d2_ij;
      }

      /*--- Calculate species laminar viscosity ---*/
      LaminarViscosity(iPoint) += (Mi/Na * gam_i) / denom;
    }
    if (ionization) {
      iSpecies = nSpecies-1;
      denom = 0.0;

      /*--- Calculate molar concentration ---*/
      Mi    = Ms[iSpecies];
      gam_i = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mi);
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        Mj    = Ms[jSpecies];
        gam_j = Primitive(iPoint,RHOS_INDEX+jSpecies) / (rho*Mj);

        /*--- Calculate "delta" quantities ---*/
        Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
            * pow(Tve, Omega11[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega11[iSpecies][jSpecies][1]*log(Tve)
            + Omega11[iSpecies][jSpecies][2]);
        d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

        /*--- Add to denominator of viscosity ---*/
        denom += gam_j*d2_ij;
      }
      LaminarViscosity(iPoint) += (Mi/Na * gam_i) / denom;
    }
  }
}

void CTNE2NSVariable ::SetThermalConductivity_GuptaYos(CConfig *config) {
  unsigned short iSpecies, jSpecies, nHeavy, nEl;
  su2double rho, T, Tve, Cvve;
  su2double *xi, *Ms, Mi, Mj, mi, mj, pi, R, RuSI, Ru, Na, kb, gam_i, gam_j, Theta_v;
  su2double denom_t, denom_r, d1_ij, d2_ij, a_ij;
  su2double ***Omega00, ***Omega11, Omega_ij;

  if (ionization) {
    cout << "SetThermalConductivity: NEEDS REVISION w/ IONIZATION" << endl;
    exit(1);
  }

  /*--- Acquire gas parameters from CConfig ---*/
  Omega00 = config->GetCollisionIntegral00();
  Omega11 = config->GetCollisionIntegral11();
  Ms      = config->GetMolar_Mass();
  xi      = config->GetRotationModes();
  if (ionization) {nHeavy = nSpecies-1;  nEl = 1;}
  else            {nHeavy = nSpecies;    nEl = 0;}

  /*--- Loop of all points ---*/
  for (unsigned long iPoint=0; iPoint<nPoint; iPoint++){

    /*--- Rename for convenience ---*/
    rho  = Primitive(iPoint,RHO_INDEX);
    T    = Primitive(iPoint,T_INDEX);
    Tve  = Primitive(iPoint,TVE_INDEX);
    Cvve = Primitive(iPoint,RHOCVVE_INDEX)/rho;
    pi   = PI_NUMBER;
    RuSI = UNIVERSAL_GAS_CONSTANT;
    Ru   = 1000.0*RuSI;
    Na   = AVOGAD_CONSTANT;
    kb   = BOLTZMANN_CONSTANT;

    /*--- Calculate mixture gas constant ---*/
    R = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      R += Ru * Primitive(iPoint,RHOS_INDEX+iSpecies)/rho;
    }

    /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
    ThermalCond    = 0.0;
    ThermalCond_ve = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

      /*--- Calculate molar concentration ---*/
      Mi      = Ms[iSpecies];
      mi      = Mi/Na;
      gam_i   = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mi);
      Theta_v = config->GetCharVibTemp(iSpecies);

      denom_t = 0.0;
      denom_r = 0.0;
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        Mj    = config->GetMolar_Mass(jSpecies);
        mj    = Mj/Na;
        gam_j = Primitive(iPoint,RHOS_INDEX+iSpecies) / (rho*Mj);

        a_ij = 1.0 + (1.0 - mi/mj)*(0.45 - 2.54*mi/mj) / ((1.0 + mi/mj)*(1.0 + mi/mj));

        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
            * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
            + Omega00[iSpecies][jSpecies][1]*log(T)
            + Omega00[iSpecies][jSpecies][2]);

        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

        /*--- Calculate the Omega^(1,1)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
            * pow(T, Omega11[iSpecies][jSpecies][0]*log(T)*log(T)
            + Omega11[iSpecies][jSpecies][1]*log(T)
            + Omega11[iSpecies][jSpecies][2]);

        /*--- Calculate "delta2_ij" ---*/
        d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

        denom_t += a_ij*gam_j*d2_ij;
        denom_r += gam_j*d1_ij;
      }

      /*--- Translational contribution to thermal conductivity ---*/
      ThermalCond(iPoint)    += (15.0/4.0)*kb*gam_i/denom_t;

      /*--- Translational contribution to thermal conductivity ---*/
      if (xi[iSpecies] != 0.0)
        ThermalCond(iPoint)  += kb*gam_i/denom_r;

      /*--- Vibrational-electronic contribution to thermal conductivity ---*/
      ThermalCond_ve(iPoint) += kb*Cvve/R*gam_i / denom_r;
    }
  }
}

void CTNE2NSVariable::SetTransportCoefficients_WBE(CConfig *config) {

  unsigned short iSpecies, jSpecies;
  su2double *Ms, Mi, Mj, M;
  su2double rho, T, Tve, RuSI, Ru, *xi;
  su2double Xs[nSpecies], conc;
  su2double Cves;
  su2double phis[nSpecies], mus[nSpecies], ks[nSpecies], kves[nSpecies];
  su2double denom, tmp1, tmp2;
  su2double **Blottner;
  su2double ***Omega00, Omega_ij;

  /*-- Loop over all points ---*/
  for (unsigned long iPoint=0; iPoint<nPoint; iPoint++){

    /*--- Rename for convenience ---*/
    rho  = Primitive(iPoint,RHO_INDEX);
    T    = Primitive(iPoint,T_INDEX);
    Tve  = Primitive(iPoint,TVE_INDEX);
    Ms   = config->GetMolar_Mass();
    xi   = config->GetRotationModes();
    RuSI = UNIVERSAL_GAS_CONSTANT;
    Ru   = 1000.0*RuSI;

    /*--- Acquire collision integral information ---*/
    Omega00 = config->GetCollisionIntegral00();

    /*--- Calculate species mole fraction ---*/
    conc = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Xs[iSpecies] = Primitive(iPoint,RHOS_INDEX+iSpecies)/Ms[iSpecies];
      conc        += Xs[iSpecies];
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Xs[iSpecies] = Xs[iSpecies]/conc;

    /*--- Calculate mixture molar mass (kg/mol) ---*/
    // Note: Species molar masses stored as kg/kmol, need 1E-3 conversion
    M = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      M += Ms[iSpecies]*Xs[iSpecies];
    M = M*1E-3;

    /*---+++                  +++---*/
    /*--- Diffusion coefficients ---*/
    /*---+++                  +++---*/

    /*--- Solve for binary diffusion coefficients ---*/
    // Note: Dij = Dji, so only loop through req'd indices
    // Note: Correlation requires kg/mol, hence 1E-3 conversion from kg/kmol
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Mi = Ms[iSpecies]*1E-3;
      for (jSpecies = iSpecies; jSpecies < nSpecies; jSpecies++) {
        Mj = Ms[jSpecies]*1E-3;

        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20/PI_NUMBER * Omega00[iSpecies][jSpecies][3]
            * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
            +  Omega00[iSpecies][jSpecies][1]*log(T)
            +  Omega00[iSpecies][jSpecies][2]);

        Dij(iPoint,iSpecies,jSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(rho*Omega_ij);
        Dij(iPoint,jSpecies,iSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(rho*Omega_ij);
      }
    }

    /*--- Calculate species-mixture diffusion coefficient --*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      denom = 0.0;
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        if (jSpecies != iSpecies) {
          denom += Xs[jSpecies]/Dij(iPoint,iSpecies,jSpecies);
        }
      }
      DiffusionCoeff(iPoint,iSpecies) = (1-Xs[iSpecies])/denom;
      //    DiffusionCoeff[iSpecies] = 0.0;
    }


    /*---+++             +++---*/
    /*--- Laminar viscosity ---*/
    /*---+++             +++---*/

    /*--- Get Blottner coefficients ---*/
    Blottner = config->GetBlottnerCoeff();

    /*--- Use Blottner's curve fits for species viscosity ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
                              Blottner[iSpecies][1])*log(T) +
          Blottner[iSpecies][2]);

    /*--- Determine species 'phi' value for Blottner model ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      phis[iSpecies] = 0.0;
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(Ms[jSpecies]/Ms[iSpecies], 0.25);
        tmp2 = sqrt(8.0*(1.0+Ms[iSpecies]/Ms[jSpecies]));
        phis[iSpecies] += Xs[jSpecies]*tmp1*tmp1/tmp2;
      }
    }

    /*--- Calculate mixture laminar viscosity ---*/
    LaminarViscosity(iPoint) = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      LaminarViscosity(iPoint) += Xs[iSpecies]*mus[iSpecies]/phis[iSpecies];

    /*---+++                +++---*/
    /*--- Thermal conductivity ---*/
    /*---+++                +++---*/

    /*--- Determine species tr & ve conductivities ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Cves = CalcCvve(Tve, config, iSpecies);
      ks[iSpecies] = mus[iSpecies]*(15.0/4.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
      kves[iSpecies] = mus[iSpecies]*Cves;
    }

    /*--- Calculate mixture tr & ve conductivities ---*/
    ThermalCond(iPoint)     = 0.0;
    ThermalCond_ve(iPoint)  = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      ThermalCond(iPoint)    += Xs[iSpecies]*ks[iSpecies]/phis[iSpecies];
      ThermalCond_ve(iPoint) += Xs[iSpecies]*kves[iSpecies]/phis[iSpecies];
    }
  }
}

bool CTNE2NSVariable::SetVorticity(void) {

  for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {

    su2double u_y = Gradient_Primitive(iPoint, VEL_INDEX, 1);
    su2double v_x = Gradient_Primitive(iPoint, VEL_INDEX+1, 0);
    su2double u_z = 0.0;
    su2double v_z = 0.0;
    su2double w_x = 0.0;
    su2double w_y = 0.0;

    if (nDim == 3) {
      u_z = Gradient_Primitive(iPoint,VEL_INDEX, 2);
      v_z = Gradient_Primitive(iPoint,VEL_INDEX+1, 2);
      w_x = Gradient_Primitive(iPoint,VEL_INDEX+2, 0);
      w_y = Gradient_Primitive(iPoint,VEL_INDEX+2, 1);
    }

    Vorticity(iPoint,0) = w_y-v_z;
    Vorticity(iPoint,1) = -(w_x-u_z);
    Vorticity(iPoint,2) = v_x-u_y;

  }
  return false;
}

bool CTNE2NSVariable::SetPrimVar_Compressible(unsigned long iPoint, CConfig *config) {

  bool nonPhys, bkup;
  unsigned short iVar;

  nonPhys = Cons2PrimVar(config, Solution[iPoint], Primitive[iPoint], dPdU[iPoint], dTdU[iPoint], dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);
  if (nonPhys) {
    for (iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);
    bkup = Cons2PrimVar(config, Solution[iPoint], Primitive[iPoint], dPdU[iPoint], dTdU[iPoint], dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);
  }

  SetVelocity2(iPoint);

  switch (config->GetKind_TransCoeffModel()) {
  case WBE:
    SetTransportCoefficients_WBE(config);
    break;
  case GUPTAYOS:
    SetDiffusionCoeff_GuptaYos(config);
    SetLaminarViscosity_GuptaYos(config);              // Requires temperature computation.
    SetThermalConductivity_GuptaYos(config);
    break;
  }

  return nonPhys;
}
