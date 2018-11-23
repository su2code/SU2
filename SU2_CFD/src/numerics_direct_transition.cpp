/*!
 * \file numerics_direct_transition.cpp
 * \brief This file contains all the convective term discretization.
 * \author A. Aranake
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/numerics_structure.hpp"
#include <limits>

CUpwSca_TransLM::CUpwSca_TransLM(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 CConfig *config)
 : CUpwScalar(val_nDim, val_nVar, config) {
}

void CUpwSca_TransLM::ExtraADPreaccIn() {

  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);
}

void CUpwSca_TransLM::FinishResidualCalc(su2double *val_residual,
                                         su2double **val_Jacobian_i,
                                         su2double **val_Jacobian_j,
                                         CConfig *config) {

  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];

  if (implicit) {
    val_Jacobian_i[0][0] = a0;    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;    val_Jacobian_i[1][1] = a0;

    val_Jacobian_j[0][0] = a1;    val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[1][0] = 0.0;    val_Jacobian_j[1][1] = a1;
  }
}

CAvgGrad_TransLM::CAvgGrad_TransLM(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   su2double *constants, bool correct_grad,
                                   CConfig *config)
 : CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) {

  sigma_intermittency = constants[6];
  sigma_Re_theta      = constants[7];
}

void CAvgGrad_TransLM::ExtraADPreaccIn() { }

void CAvgGrad_TransLM::FinishResidualCalc(su2double *val_residual,
                                          su2double **Jacobian_i,
                                          su2double **Jacobian_j,
                                          CConfig   *config) {

  /*--- Compute mean effective viscosity for both equations. ---*/
  const su2double diff_i_int   =  Laminar_Viscosity_i + Eddy_Viscosity_i /sigma_intermittency;
  const su2double diff_j_int   =  Laminar_Viscosity_j + Eddy_Viscosity_j /sigma_intermittency;
  const su2double diff_i_theta = (Laminar_Viscosity_i + Eddy_Viscosity_i)*sigma_Re_theta;
  const su2double diff_j_theta = (Laminar_Viscosity_j + Eddy_Viscosity_j)*sigma_Re_theta;

  const su2double diff_int   = 0.5*(diff_i_int   + diff_j_int);   // Could instead use weighted average!
  const su2double diff_theta = 0.5*(diff_i_theta + diff_j_theta);

  /*--- Compute the values of the residuals on this edge. */
  val_residual[0] = diff_int  *Proj_Mean_GradTurbVar[0];
  val_residual[1] = diff_theta*Proj_Mean_GradTurbVar[1];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_int*proj_vector_ij/Density_i;
    Jacobian_i[0][1] =  0.0;
    Jacobian_i[1][0] =  0.0;
    Jacobian_i[1][1] = -diff_theta*proj_vector_ij/Density_i;

    Jacobian_j[0][0] = diff_int*proj_vector_ij/Density_j;
    Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;
    Jacobian_j[1][1] = diff_theta*proj_vector_ij/Density_j;
  }
}

CSourcePieceWise_TransLM::CSourcePieceWise_TransLM(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   su2double* constants,
                                                   CConfig *config)
 : CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Store the closure constants of the LM model in a more readable way. ---*/
  ca1          = constants[0];
  ca2          = constants[1];
  ce1          = constants[2];
  ce2          = constants[3];
  cthetat      = constants[4];
  Flength_CF   = constants[8];
  C_Fonset1_CF = constants[9];
  CHe_max      = constants[10];
}

CSourcePieceWise_TransLM::~CSourcePieceWise_TransLM(void) { }

su2double CSourcePieceWise_TransLM::GetREth(const su2double var_tu) {

  const su2double tu = 100.0*var_tu; // Turbulence intensity in percents.

  if (tu <= 1.3)
    return (1173.9967604078363 - 589.428*tu + 0.2196/(tu*tu)); // Modified the original value of 1173.51 to get better continuity at tu=1.3.
  else
    return (331.5*pow(tu-0.5658, -0.671));
}

su2double CSourcePieceWise_TransLM::GetREth_crit(const su2double var_Re_theta) {

  su2double Recrit;

  /* Menter correlation, also mentioned at https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html.
     The coefficients are slightly altered to have a better continuity at the interface. */
  if(var_Re_theta <= 1870.0) {

    /* Quartic expression is used. Define the constants. */
    const su2double a0 = -3.96035;
    const su2double a1 =  1.0120656;
    const su2double a2 = -8.68230e-4;
    const su2double a3 =  6.96506e-7;
    const su2double a4 = -1.74105e-10;

    /* Compute the value of Recrit. */
    const su2double val1 = var_Re_theta;
    const su2double val2 = val1*val1;
    const su2double val3 = val1*val2;
    const su2double val4 = val2*val2;

    Recrit = a0 + a1*val1 + a2*val2 + a3*val3 + a4*val4;

  } else {

    /* Use correlation valid for Re_theta larger than 1870.0. */
    Recrit = var_Re_theta - (591.926884931 + 0.482*(var_Re_theta-1870.0));
  }

  /* Malan correlation. */
  // Recrit = min(0.615*var_Re_theta+61.5, var_Re_theta);

  /* Return the value of Recrit. */
  return Recrit;
}

su2double CSourcePieceWise_TransLM::GetFlength(const su2double var_Re_theta) {

  su2double Flength;

  /* Menter correlation, also mentioned at https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html.
     The coefficients are slightly changed to have a better continuity at the interfaces. */
  if(var_Re_theta < 400.0) {
    Flength = 39.825535395 - 1.1927e-2*var_Re_theta - 1.32567e-4*var_Re_theta*var_Re_theta;
  }
  else if(var_Re_theta < 596.0) {
    Flength = 263.408015395 - 1.23939*var_Re_theta + 1.94548e-3*var_Re_theta*var_Re_theta
            - 1.01695e-6*var_Re_theta*var_Re_theta*var_Re_theta;
  }
  else if(var_Re_theta < 1200.0) {
    Flength = 0.5 - 3.0e-4*(var_Re_theta - 596.0);
  }
  else {
    Flength = 0.3188;
  }
  
  /* Malan correlation. */
  /* const su2double tmp1 = 7.168 - 0.01173*var_Re_theta;
     const su2double tmp2 = exp(tmp1) + 0.5;
     Flength = min(tmp2, 300.0);
  */

  /* Return the value of Flength. */
  return Flength;
}

void CSourcePieceWise_TransLM::ComputeResidual(su2double *val_residual,
                                               su2double **val_Jacobian_i,
                                               su2double **val_Jacobian_j,
                                               CConfig   *config) {

  /*--- Determine whether a Spalart-Allmaras type turbulence model is used. ---*/
  const unsigned short turbModel = config->GetKind_Turb_Model();
  const bool SA_turb  = (turbModel == SA) || (turbModel == SA_NEG) || (turbModel == SA_E) ||
                        (turbModel == SA_COMP) || (turbModel == SA_E_COMP);

  /*--- Minimum and maximum value of lamt. The recommended values are -0.1 and
        0.1, see https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html. ---*/
  //const su2double lamt_min = -0.15;
  //const su2double lamt_max =  0.25;
  const su2double lamt_min = -0.1;
  const su2double lamt_max =  0.1;

  /*--- Minimum and maximum value of the turbulence intensity. The recommended
        value for the minimum is 0.00027, see
        see https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html. ---*/
  const su2double turbIntensity_min = 0.00027;
  const su2double turbIntensity_max = 1.0;

  /*--- Maximum number of iterations in the Newton algorithm for Re_theta_eq. ---*/
  const unsigned short maxIt = 25;

  /*--- Get the infinity state, which is needed for some of the
        cross flow instability models. ---*/
  const su2double Density_Inf   = config->GetDensity_FreeStreamND();
  const su2double Pressure_Inf  = config->GetPressure_FreeStreamND();
  const su2double *Velocity_Inf = config->GetVelocity_FreeStreamND();

  /*--- Get the specific heat ratio. ---*/
  const su2double Gamma = config->GetGamma();
  const su2double Gamma_Minus_One = Gamma - 1.0;

  /*--- More readable storage of some variables. ---*/
  const su2double rho    = V_i[nDim+2];
  const su2double p      = V_i[nDim+1];
  const su2double *vel   = V_i + 1;
  const su2double muLam  = Laminar_Viscosity_i;
  const su2double muTurb = Eddy_Viscosity_i;
  const su2double S      = StrainMag_i;
  const su2double Vort   = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1]
                         +      Vorticity_i[2]*Vorticity_i[2]);
  const su2double dist   = dist_i;

  const su2double kine  = SA_turb ? 0.0 : TurbVar_i[0];
  const su2double omega = SA_turb ? 1.0 : TurbVar_i[1];

  const su2double intermittency = TransVar_i[0];
  const su2double Re_theta      = TransVar_i[1];

  /*--- Compute the velocity magnitude, velocity magnitude squared
        of the free-stream and helicity. ---*/
  su2double vel2 = 0.0, velInf2 = 0.0, helicity = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; ++iDim) {
    vel2     += vel[iDim]*vel[iDim];
    helicity += vel[iDim]*Vorticity_i[iDim];
    velInf2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  vel2 = max(vel2, 1.e-25);
  const su2double velMag = sqrt(vel2);

  /*--- Compute the turbulence intensity. For Spalart-Allmaras type turbulence
        models the turbulence intensity of the free-stream is taken. For numerical
        robustness, the turbulence intensity is clipped within practical limits. ---*/
  su2double turbIntensity         = sqrt(2.0*kine/(3.0*vel2));
  if(dist < 1.e-10) turbIntensity = 0.0;
  if( SA_turb ) turbIntensity     = config->GetTurbulenceIntensity_FreeStream();

  turbIntensity = min(max(turbIntensity, turbIntensity_min), turbIntensity_max);

  /*--- Compute the value of dUds, which is the acceleration of the
        total velocity along streamwise direction. ---*/
  su2double dUds = 0.0;
  for(unsigned short iDim=0; iDim<nDim; ++iDim) {
    for(unsigned short jDim=0; jDim<nDim; ++jDim)
      dUds += vel[iDim]*vel[jDim]*PrimVar_Grad_i[iDim+1][jDim];
  }
  dUds /= vel2;

  /*--- Computation of the equilibrium value of Re_theta. As the equation
        for Re_theta_eq is an implicit equation, an iterative algorithm is
        needed to compute it. ---*/
  const su2double Re_theta_far = GetREth(turbIntensity);
  const su2double termLam      = muLam*dUds/(rho*vel2);

  su2double Re_theta_eq = Re_theta_far; /* Initial guess. */

  if(dUds <= 0.0) {

    /*--- Negative value of dUds implies a negative value of lambda_theta
          and a different correlation must be taken than for positive values,
          see https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html. ---*/
    const su2double val    = 100.0*turbIntensity/1.5;
    const su2double TITerm = exp(-val*sqrt(val));

    /* Constants in the correlation for lambda_theta. */
    const su2double a1 =  12.986*TITerm;
    const su2double a2 = 123.66 *TITerm;
    const su2double a3 = 405.689*TITerm;

    /* Determine the lower limit of Re_theta_eq, which corresponds
       lambda_theta = lamt_min. */
    su2double lamt  = lamt_min;
    su2double lamt2 = lamt*lamt;
    su2double lamt3 = lamt*lamt2;
    su2double F     = 1.0 + a1*lamt + a2*lamt2 + a3*lamt3;
    su2double dF;

    const su2double Re_theta_low = Re_theta_far*F;

    /*--- Newton algorithm to compute Re_theta_eq. ---*/
    unsigned short iter;
    su2double deltaReOld = 1.0;
    for(iter=0; iter<maxIt; ++iter) {

      /* Value of lamt. For stability this term is clipped. */
      lamt = Re_theta_eq*Re_theta_eq*termLam;
      if(lamt < lamt_min) {
        lamt        = lamt_min;
        Re_theta_eq = Re_theta_low;
      }

      /* Compute the value of the function F(lamt) and its derivative
         w.r.t. lamt. */
      lamt2 = lamt*lamt;
      lamt3 = lamt*lamt2;
      F     = 1.0 + a1*lamt + a2*lamt2 + a3*lamt3;
      dF    = a1 + 2.0*a2*lamt + 3.0*a3*lamt2;

      /* Compute the function for which the root must be determined
         as well as its derivative w.r.t. Re_theta_eq. */
      const su2double G  = Re_theta_eq - Re_theta_far*F;
      const su2double dG = 1.0 - Re_theta_far*dF*2.0*Re_theta_eq*termLam;

      /* Update the value of Re_theta_eq. Store the old value, because
         the new value may need to be clipped. */
      const su2double Re_theta_eq_Old = Re_theta_eq;
      Re_theta_eq -= G/dG;

      /* Clip Re_theta_eq for stability and determine the actual change. */
      Re_theta_eq = max(min(Re_theta_eq, Re_theta_far), Re_theta_low);

      const su2double deltaRe = fabs(Re_theta_eq - Re_theta_eq_Old);

      /* Exit criterion, which takes finite precision into account. */
      if((deltaRe <= 1.e-2) && (deltaRe >= deltaReOld)) break;
      deltaReOld = deltaRe;
    }

    /* Terminate if the Newton algorithm did not converge. */
    if(iter == maxIt)
      SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);
  }
  else {

    /*--- Positive value of dUds implies a positive value of lambda_theta
          and a different correlation must be taken than for negative values,
          see https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html. ---*/
    const su2double val    = 100.0*turbIntensity/0.5;
    const su2double TITerm = exp(-val);

    /* Constants in the correlation for lambda_theta. */
    const su2double a1 =   0.275*TITerm;
    const su2double a2 = -35.0;

    /* Determine the upper limit of Re_theta_eq, which corresponds
       lambda_theta = lamt_max. */
    su2double lamt = lamt_max;
    su2double valExp = exp(a2*lamt);
    su2double F      = 1.0 + a1*(1.0 - valExp);
    su2double dF;

    const su2double Re_theta_upp = Re_theta_far*F;

    /*--- Newton algorithm to compute Re_theta_eq. ---*/
    unsigned short iter;
    su2double deltaReOld = 1.0;
    for(iter=0; iter<maxIt; ++iter) {

      /* Value of lamt. For stability this term is clipped. */
      lamt = Re_theta_eq*Re_theta_eq*termLam;
      if(lamt > lamt_max) {
        lamt = lamt_max;
        Re_theta_eq = Re_theta_upp;
      }

      /* Compute the value of the function F(lamt) and its derivative
         w.r.t. lamt. */
      valExp = exp(a2*lamt);
      F      = 1.0 + a1*(1.0 - valExp);
      dF     = -a1*a2*valExp;

      /* Compute the function for which the root must be determined
         as well as its derivative w.r.t. Re_theta_eq. */
      const su2double G  = Re_theta_eq - Re_theta_far*F;
      const su2double dG = 1.0 - Re_theta_far*dF*2.0*Re_theta_eq*termLam;

      /* Update the value of Re_theta_eq. Store the old value, because
         the new value may need to be clipped. */
      const su2double Re_theta_eq_Old = Re_theta_eq;
      Re_theta_eq -= G/dG;

      /* Clip Re_theta_eq for stability and determine the actual change. */
      Re_theta_eq = max(min(Re_theta_eq, Re_theta_upp), Re_theta_far);

      const su2double deltaRe = fabs(Re_theta_eq - Re_theta_eq_Old);

      /* Exit criterion, which takes finite precision into account. */
      if((deltaRe <= 1.e-2) && (deltaRe >= deltaReOld)) break;
      deltaReOld = deltaRe;
    }

    /* Terminate if the Newton algorithm did not converge. */
    if(iter == maxIt)
      SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);
  }

  /*--- Compute the other Reynolds numbers that appear in the formulation. ---*/
  const su2double Rev           = rho*S*dist*dist/muLam;
  const su2double Re_theta_crit = GetREth_crit(Re_theta);
  const su2double Re_omega      = rho*omega*dist*dist/muLam;

  /*--- Compute the value of Fonset. ---*/
  const su2double RT = SA_turb ? muTurb/muLam : rho*kine/(muLam*omega);

  const su2double Fonset1 = Rev/(2.193*Re_theta_crit);
  su2double       tmp     = Fonset1*Fonset1*Fonset1*Fonset1;
  const su2double Fonset2 = min(max(Fonset1,tmp), 2.0);
  tmp                     = 1.0 - 0.064*RT*RT*RT;          // 2.5^(-3) = 0.064.
  const su2double Fonset3 = max(tmp, 0.0);
  const su2double Fonset  = max(Fonset2-Fonset3, 0.0);

  /*--- Compute the value of Fturb. ---*/
  const su2double Fturb = exp(-0.00390625*RT*RT*RT*RT);       // 4^(-4) = 0.00390625.

  /*--- Compute the values of Fsublayer and Fwake. When a
        Spalart-Allmaras type model is used, these value
        are set to be zero and one, respectively. ---*/
  tmp                       = 2.5e-5*Re_omega*Re_omega;
  const su2double Fsublayer = SA_turb ? 0.0 : exp(-tmp);
  tmp                       = 1.e-10*Re_omega*Re_omega;
  const su2double Fwake     = SA_turb ? 1.0 : exp(-tmp);

  /*--- Compute the empirical correlation for Flength, which
        controls the length of the transition region. It is,
        possibly, corrected by Fsublayer. ---*/
  const su2double Flength = (1.0-Fsublayer)*GetFlength(Re_theta)
                          + 40.0*Fsublayer;

  /*--- Compute the value of Ftheta_t. */
  tmp                  = rho*vel2/(375.0*Vort*muLam*Re_theta);
  const su2double val1 = Fwake*exp(-tmp*tmp*tmp*tmp);
  const su2double val2 = (ce2*intermittency-1.0)/(ce2-1.0);
  const su2double val3 = 1.0 - val2*val2;

  const su2double Ftheta_t = min(max(val1,val3),1.0);

  /*--- Compute the helicity Reynolds number. ---*/
  const su2double Re_helicity = rho*dist*dist*helicity/(muLam*velMag);

  /*--- Compute the velocity at the edge of the boundary layer. Here it is
        assumed that the isentropic relation for total pressure is valid in
        the outer flow and the pressure is constant through the boundary
        layer. Note that this relation is not valid when motion is present. ---*/
  const su2double pRatio = p/Pressure_Inf;
  tmp  = 1.0 - pow(pRatio,(Gamma_Minus_One/Gamma));
  tmp *= 2.0*Gamma*Pressure_Inf/(Gamma_Minus_One*Density_Inf);
  tmp += velInf2;

  su2double uEdge = sqrt(max(tmp,0.0));
  uEdge = max(uEdge, 0.01*velMag);

  /*--- Compute the dot product of the velocity and the pressure gradient.
        Store this in duEdgeDs, because later on this value is transformed
        to the gradient of uEdge in the direction of the streamline. ---*/
  su2double duEdgeDs = 0.0;
  for(unsigned short iDim=0; iDim<nDim; ++iDim)
    duEdgeDs += vel[iDim]*PrimVar_Grad_i[nDim+1][iDim];

  /*--- Multiply the dot product of velocity and pressure gradient with
        the appropriate variables to obtain the gradient of the edge
        velocity in the direction of the streamline. ---*/
  tmp = pow(pRatio, 1.0/Gamma);
  duEdgeDs *= -1.0/(tmp*Density_Inf*uEdge*velMag);

  /*--- Compute the length scale that represents the momentum thickness
        where the helicity Reynolds number reaches its maximum value
        in the boundary layer and the modified pressure gradient parameter
        that is needed to compute the shape factor H12. ---*/
  const su2double lenHelMax = 2.0*dist/(15.0*CHe_max);
  const su2double lamPlus   = rho*lenHelMax*lenHelMax*duEdgeDs/muLam;

  /*--- Compute the empirical correlation for the helicity Reynolds number
        at transition onset. ---*/
  const su2double lamPlus2 = lamPlus*lamPlus;
  const su2double lamPlus3 = lamPlus*lamPlus2;
  const su2double lamPlus4 = lamPlus*lamPlus3;

  tmp = -8838.4*lamPlus4 + 1105.1*lamPlus3 - 67.962*lamPlus2 + 17.574*lamPlus + 2.0593;
  const su2double H12 = 4.02923 - sqrt(max(tmp,0.0));

  tmp = -456.83*H12 + 1332.7;
  const su2double Re_helicity_onset = max(tmp, 150.0);

  /*--- Compute the value of Fonset_CF. ---*/
  const su2double Fonset1_CF = Re_helicity/(C_Fonset1_CF*Re_helicity_onset);

  tmp = Fonset1_CF*Fonset1_CF*Fonset1_CF*Fonset1_CF;
  const su2double Fonset2_CF = min(max(Fonset1_CF,tmp), 2.0);

  tmp = 1.0 - 0.125*RT*RT*RT;
  const su2double Fonset3_CF = max(tmp, 0.0);

  tmp = Fonset2_CF - Fonset3_CF;
  const su2double Fonset_CF = max(tmp, 0.0);

  /*--- Compute the production and destruction term of the intermittency
        equation. Note that Flength_CF is set to zero when the crossflow
       	instability part is switched off. ---*/
  const su2double termPInter     = Flength*sqrt(intermittency*Fonset)
                                 + Flength_CF*sqrt(intermittency*Fonset_CF);
  const su2double Pintermittency = termPInter*ca1*rho*S*(1.0 - ce1*intermittency);
  const su2double Eintermittency = ca2*rho*Vort*Fturb
                                 * intermittency*(ce2*intermittency - 1.0);

  /*--- Compute the production/destruction term of the Re_theta equation. ---*/
  const su2double PRe_theta = cthetat*rho*rho*vel2*(Re_theta_eq - Re_theta)
                            * (1.0 - Ftheta_t)/(500.0*muLam);

  /*--- Compute the source terms and the approximate Jacobians for
        both equations. ---*/
  val_residual[0] = Volume*(Pintermittency - Eintermittency);
  val_residual[1] = Volume*PRe_theta;

  val_Jacobian_i[0][0] = -Volume*(2.0*ca2*ce2*Vort*intermittency*Fturb
                       +          ca1*ce1*S*termPInter);
  val_Jacobian_i[0][1] =  0.0;
  val_Jacobian_i[1][0] =  0.0;
  val_Jacobian_i[1][1] = -Volume*cthetat*rho*vel2/(500.0*muLam);
}
