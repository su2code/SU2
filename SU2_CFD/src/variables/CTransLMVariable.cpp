/*!
 * \file CTransLMVariable.cpp
 * \brief Definition of the solution fields.
 * \author A. Aranake, A. Rausa, M. Cerabona
 * \version 7.4.0 "Blackbird"
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

#include "../../include/variables/CTransLMVariable.hpp"
#include <cmath>


CTransLMVariable::CTransLMVariable(su2double intermittency, su2double Re_theta, const su2double* constants, unsigned long npoint, unsigned long ndim,
                                         unsigned long nvar, CConfig *config) :
                                         CTurbVariable(npoint, ndim, nvar, config) {

    for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {
        Solution(iPoint,0) = intermittency;
        Solution(iPoint,1) = Re_theta;
    }

    Solution_Old = Solution;  


    F_onset.resize(nPoint) = su2double(0.0);
    F_turb.resize(nPoint) = su2double(0.0);
    FReattach.resize(nPoint) = su2double(0.0);
    F_length1.resize(nPoint) = su2double(0.0);
    F_length.resize(nPoint) = su2double(0.0);
    F_sublayer.resize(nPoint) = su2double(0.0);
    F_thetat.resize(nPoint) = su2double(0.0);
    F_wake.resize(nPoint) = su2double(0.0);
    F_lambda.resize(nPoint) = su2double(0.0);


    reV.resize(nPoint) = su2double(0.0); 
    R_T.resize(nPoint) = su2double(0.0);
    rew.resize(nPoint) = su2double(0.0);
    rethetac.resize(nPoint) = su2double(0.0);
    rethetat_eq.resize(nPoint) = su2double(0.0);

    T.resize(nPoint) = su2double(0.0);
    delta_param.resize(nPoint) = su2double(0.0);
    lambda_theta.resize(nPoint) = su2double(0.0);
    Turb_Intens.resize(nPoint) = su2double(0.0);
    dU_ds.resize(nPoint) = su2double(0.0);
    thetat.resize(nPoint) = su2double(0.0);
    gamma_sep.resize(nPoint) = su2double(0.0);
    gamma_eff.resize(nPoint) = su2double(0.0);


    // Varibles needed in case SA has to be converted to SST
    k.resize(nPoint) = su2double(0.0);
    w.resize(nPoint) = su2double(0.0);


    // Variables needed for LM2015 model
    F_thetat2.resize(nPoint) = su2double(0.0);
    ReThetat_SCF.resize(nPoint) = su2double(0.0);
    thetat_SCF.resize(nPoint) = su2double(0.0);


    TurbSA = TurbModelFamily(config->GetKind_Turb_Model()) == TURB_FAMILY::SA;
    TurbKW = TurbModelFamily(config->GetKind_Turb_Model()) == TURB_FAMILY::KW;
    TransModel = config->GetKind_Trans_Model();

    if(TurbSA)

    if(TurbSA)
      Turb_Intens = config->GetTurbulenceIntensity_FreeStream()*100;

    Trans_Correlation = config->GetKind_Trans_Correlation();

    if(Trans_Correlation == TURB_TRANS_CORRELATION::DEFAULT && TurbSA)
      Trans_Correlation = TURB_TRANS_CORRELATION::MALAN;
    if(Trans_Correlation == TURB_TRANS_CORRELATION::DEFAULT && TurbKW)
      Trans_Correlation = TURB_TRANS_CORRELATION::MENTER_LANGTRY;


}


void CTransLMVariable::Set_kAndw(unsigned long iPoint, su2double VorticityMag, su2double val_velocity, su2double val_eddy_viscosity, su2double StrainMag, su2double val_dist,
               su2double val_viscosity, CConfig *config){

  if(TurbSA && config->GetConvertSA2SST()) {
    su2double BetaStar = 0.09;
    w(iPoint) = VorticityMag / sqrt(BetaStar);

    su2double FirstTerm = 2.0 * sqrt(val_eddy_viscosity / w(iPoint)) / (BetaStar * val_dist);
    su2double SecondTerm = 500 * val_viscosity / (val_dist * val_dist * w(iPoint));
    su2double arg2 = max(FirstTerm, SecondTerm);
    su2double F2 = tanh(arg2 * arg2);

    su2double a1 = 0.31;
    k(iPoint) = val_eddy_viscosity * max(w(iPoint), VorticityMag * F2 / a1);


//    Turb_Intens(iPoint) = 100*sqrt(2*k(iPoint)/3)/val_velocity;
//    Turb_Intens(iPoint) = max(Turb_Intens(iPoint), 0.027);

  }




}

void CTransLMVariable::ReThetaC_Correlations(unsigned long iPoint){

  su2double ReTheta_Tilde = Solution(iPoint, 1);

  switch (Trans_Correlation) {
    case TURB_TRANS_CORRELATION::MALAN: {
      rethetac(iPoint) = min(0.615 * ReTheta_Tilde + 61.5, ReTheta_Tilde);
      break;
    }

    case TURB_TRANS_CORRELATION::SULUKSNA: {
      rethetac(iPoint) = min(0.1 * exp(-0.0022 * ReTheta_Tilde + 12), 300.0);
      break;
    }

    case TURB_TRANS_CORRELATION::KRAUSE: {
      rethetac(iPoint) = 0.91 * ReTheta_Tilde + 5.32;
      break;
    }

    case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
      su2double FirstTerm = 4.45 * pow(Turb_Intens(iPoint), 3);
      su2double SecondTerm = 5.7 * pow(Turb_Intens(iPoint), 2);
      su2double ThirdTerm = 1.37 * pow(Turb_Intens(iPoint), 1);
      rethetac(iPoint) = (FirstTerm - SecondTerm + ThirdTerm + 0.585) * ReTheta_Tilde;
      break;
    }

    case TURB_TRANS_CORRELATION::MEDIDA: {
      rethetac(iPoint) = 0.62 * ReTheta_Tilde;
      break;
    }

    case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {

      if (ReTheta_Tilde <= 1870) {
        su2double FirstTerm = (-396.035 * pow(10, -2));
        su2double SecondTerm = (10120.656 * pow(10, -4)) * ReTheta_Tilde;
        su2double ThirdTerm = (-868.230 * pow(10, -6)) * pow(ReTheta_Tilde, 2);
        su2double ForthTerm = (696.506 * pow(10, -9)) * pow(ReTheta_Tilde, 3);
        su2double FifthTerm = (-174.105 * pow(10, -12)) * pow(ReTheta_Tilde, 4);
        rethetac(iPoint) = FirstTerm + SecondTerm + ThirdTerm + ForthTerm + FifthTerm;
      } else {
        rethetac(iPoint) = ReTheta_Tilde - (593.11 + 0.482 * (ReTheta_Tilde - 1870.0));
      }

      break;
    }

  }

  rethetac(iPoint) = max(rethetac(iPoint), 1e-20);

}

void CTransLMVariable::FLength_Correlations(unsigned long iPoint){

  su2double ReTheta_Tilde = Solution(iPoint, 1);

  switch (Trans_Correlation) {
    case TURB_TRANS_CORRELATION::MALAN: {
      F_length1(iPoint) = min(exp(7.168 - 0.01173 * ReTheta_Tilde) + 0.5, 300.0);
      break;
    }

    case TURB_TRANS_CORRELATION::SULUKSNA: {
      su2double FirstTerm = -pow(0.025 * ReTheta_Tilde, 2) + 1.47 * ReTheta_Tilde - 120.0;
      F_length1(iPoint) = min(max(FirstTerm, 125.0), ReTheta_Tilde);
      break;
    }

    case TURB_TRANS_CORRELATION::KRAUSE: {
      F_length1(iPoint) = 3.39 * ReTheta_Tilde + 55.03;
      break;
    }

    case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
      su2double FirstTerm = 0.171 * pow(Turb_Intens(iPoint), 2);
      su2double SecondTerm = 0.0083 * pow(Turb_Intens(iPoint), 1);
      F_length1(iPoint) = (FirstTerm - SecondTerm + 0.0306);
      break;
    }

    case TURB_TRANS_CORRELATION::MEDIDA: {
      F_length1(iPoint) = 40;
      break;
    }

    case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {
      if (ReTheta_Tilde < 400) {
        F_length1(iPoint) = 39.8189 + (-119.270 * pow(10, -4)) * ReTheta_Tilde +
                            (-132.567 * pow(10, -6)) * ReTheta_Tilde * ReTheta_Tilde;
      } else if (ReTheta_Tilde < 596) {
        F_length1(iPoint) = 263.404 + (-123.939 * pow(10, -2)) * ReTheta_Tilde +
                            (194.548 * pow(10, -5)) * pow(ReTheta_Tilde, 2) +
                            (-101.695 * pow(10, -8)) * pow(ReTheta_Tilde, 3);
      } else if (ReTheta_Tilde < 1200) {
        F_length1(iPoint) = 0.5 - (3.0 * pow(10, -4)) * (ReTheta_Tilde - 596.0);
      } else {
        F_length1(iPoint) = 0.3188;
      }
      break;
    }

  }

}


void CTransLMVariable::SetF_onset(unsigned long iPoint, su2double val_density, su2double StrainMag, su2double val_dist,
                                  su2double val_viscosity, su2double *TurbVars, su2double val_eddy_viscosity, su2double VorticityMag, su2double val_velocity, CConfig *config) {


  reV(iPoint) = val_density*StrainMag*val_dist*val_dist/val_viscosity;

  // If model is coupled with SST then TI can be locally computed
  if (TurbKW){
    Turb_Intens(iPoint) = 100*sqrt(2*TurbVars[0]/3)/val_velocity;
    Turb_Intens(iPoint) = max(Turb_Intens(iPoint), 0.027);
  }



  // Computing ReTheta_C from correlations
  // For LM-SST couple, the preferred one is LANGTRY_MENTER
  // For LM-SA couple, the preferred one is MALAN
  // Change these in the configuration file!
  ReThetaC_Correlations(iPoint);



  //    cout << "Dopo rethetac" << endl;

  /*--- F_onset ---*/
  su2double F_onset1 = reV(iPoint)/(2.193*rethetac(iPoint));

  // If model is coupled with SA then R_T is expressed differently
  if (TurbSA)
    R_T(iPoint) = val_eddy_viscosity / val_viscosity;

  // Original formulation
  if (TurbKW)
    R_T(iPoint) = val_density*TurbVars[0]/(val_viscosity*TurbVars[1]);


  if (TurbSA && config->GetConvertSA2SST())
    R_T(iPoint) = val_density*k(iPoint)/(val_viscosity*w(iPoint));

//  cout << R_T(iPoint) << endl;

  su2double F_onset2 = 0.0;
  su2double F_onset3 = 0.0;

  if (TurbSA) {
    F_onset2 = min(max(F_onset1, pow(F_onset1, 4)), 4.0);
    F_onset3 = max(2.0 - pow(0.4 * R_T(iPoint), 3), 0.0);

  }
  if (TurbKW) {
    F_onset2 = min(max(F_onset1, pow(F_onset1, 4)), 2.0);
    F_onset3 = max(1.0 - pow(0.4 * R_T(iPoint), 3), 0.0);
  }

  F_onset(iPoint) = max(F_onset2-F_onset3,0.0);

}

void CTransLMVariable::SetF_length(su2double *Velocity, unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_dist,
                                   su2double val_viscosity, CConfig *config) {

  // Computing F_Length from correlations
  // For LM-SST couple, the preferred one is LANGTRY_MENTER
  // For LM-SA couple, the preferred one is MALAN
  // Change these in the configuration file!
  FLength_Correlations(iPoint);

  if(TurbSA){

    F_length(iPoint) = F_length1(iPoint);
  }


  if (TurbKW) {

    rew(iPoint) =
        val_density * TurbVars[1] * val_dist * val_dist / val_viscosity;

    F_sublayer(iPoint) = exp(-pow(0.005 * rew(iPoint), 2));

    F_length(iPoint) = F_length1(iPoint) * (1 - F_sublayer(iPoint)) + 40.0 * F_sublayer(iPoint);

  }

  if(TurbSA && config->GetConvertSA2SST()) {

    F_length(iPoint) = F_length1(iPoint);

  }

}

void CTransLMVariable::SetF_turb(unsigned long iPoint) {

  F_turb(iPoint) = exp(-pow(0.25*R_T(iPoint),4));

}

void CTransLMVariable::Setrethetat_eq(unsigned long iPoint, su2double *Velocity, su2double val_velocity, CMatrixView<const su2double> Velocity_Gradient,
                                      su2double *TurbVars, su2double val_viscosity, su2double val_density, CConfig *config) {


  // Compute derivatives of velocity field
  su2double du_dx = 0.0;
  for(auto iDim = 0u; iDim < nDim; iDim++){
    du_dx += 2 * Velocity[iDim] * Velocity_Gradient[iDim][0];
  }
  du_dx = du_dx * 0.5 / val_velocity;

  su2double du_dy = 0.0;
  for(auto iDim = 0u; iDim < nDim; iDim++){
    du_dy += 2 * Velocity[iDim] * Velocity_Gradient[iDim][1];
  }
  du_dy = du_dy * 0.5 / val_velocity;

  su2double du_dz = 0.0;
  if(nDim == 3) {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      du_dz += 2 * Velocity[iDim] * Velocity_Gradient[iDim][2];
    }
    du_dz = du_dz * 0.5 / val_velocity;
  }



  dU_ds(iPoint) = Velocity[0] * du_dx / val_velocity + Velocity[1] * du_dy / val_velocity;
  if (nDim==3) {
    dU_ds(iPoint) += Velocity[2] * du_dz / val_velocity;
  }


  // Initial guess is made with lambda = 0
  F_lambda(iPoint) = 1.0;
  su2double toll = 1e-5;
  su2double error = toll+1.0;
  su2double rethetat_eq_old = 20.0;   
  int nMax = 100;


  int iter;
  for (iter=0;iter<nMax && error > toll;iter++) {  

    thetat(iPoint) = rethetat_eq_old*val_viscosity/(val_density*val_velocity);
    lambda_theta(iPoint) = val_density*thetat(iPoint)*thetat(iPoint)*dU_ds(iPoint)/val_viscosity;
    
    // Limit lambda_theta between -0.1 and 0.1
    lambda_theta(iPoint) = max(-0.1, lambda_theta(iPoint));
    lambda_theta(iPoint) = min(0.1, lambda_theta(iPoint));

    if (lambda_theta(iPoint)<=0.0) {
      su2double FirstTerm = 12.986*lambda_theta(iPoint);
      su2double SecondTerm = 123.66*pow(lambda_theta(iPoint),2);
      su2double ThirdTerm = 405.689*pow(lambda_theta(iPoint),3);
      F_lambda(iPoint) = 1 + ( FirstTerm + SecondTerm + ThirdTerm ) * exp(- pow(Turb_Intens(iPoint)/1.5,1.5));
    }
    else {
      su2double FirstTerm = exp(-35.0*lambda_theta(iPoint));
      F_lambda(iPoint) = 1.0 + 0.275*(1.0 - FirstTerm) * exp(-2.0*Turb_Intens(iPoint));
    }

    if (Turb_Intens(iPoint)<=1.3) {
      su2double FirstTerm = 589.428*Turb_Intens(iPoint);
      su2double SecondTerm = 0.2196/(Turb_Intens(iPoint)*Turb_Intens(iPoint));
      rethetat_eq(iPoint) = (1173.51 - FirstTerm + SecondTerm) * F_lambda(iPoint);
    }
    else {
      rethetat_eq(iPoint) = 331.5*pow(Turb_Intens(iPoint)-0.5658,-0.671)*F_lambda(iPoint);
    }

    rethetat_eq(iPoint) = max(20.0, rethetat_eq(iPoint));   

    error = abs(rethetat_eq(iPoint) - rethetat_eq_old)/rethetat_eq_old;

    rethetat_eq_old = rethetat_eq(iPoint);


  }



}

void CTransLMVariable::SetT(unsigned long iPoint, su2double val_velocity, su2double val_viscosity, su2double val_density, su2double val_localGridLength, su2double val_mu_T) {


  T(iPoint) = 500*val_viscosity/(val_density*val_velocity*val_velocity);   

  if(TransModel == TURB_TRANS_MODEL::LM2015){
    T(iPoint) = min(T(iPoint), val_density*val_localGridLength*val_localGridLength/(val_viscosity + val_mu_T));
  }

  T(iPoint) = max(T(iPoint), 1e-20);

}

void CTransLMVariable::SetF_thetat(unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_dist,
                                   su2double val_viscosity, su2double val_vort, su2double val_velocity, su2double* constants, CConfig *config) {


  su2double theta_BL = Solution(iPoint, 1) * val_viscosity / (val_density * val_velocity);
  su2double delta_BL = 15.0 * theta_BL / 2.0;
  delta_param(iPoint) = 50.0 * val_vort * val_dist * delta_BL / val_velocity;
  delta_param(iPoint) = max(delta_param(iPoint), 1e-20);


  // If model is coupled with SA then F_wake = 1
  if (TurbSA)
    F_wake(iPoint) = 1.0;

  if (TurbKW)
    F_wake(iPoint) = exp(-pow(0.00001*rew(iPoint),2));

  if (TurbSA && config->GetConvertSA2SST())
    F_wake(iPoint) = 1.0;


  su2double c_e2 = constants[3];
  su2double FirstMaxTerm = F_wake(iPoint) * exp(-pow(val_dist/delta_param(iPoint), 4));
  su2double SecondMaxTerm = 1 - pow((c_e2 * Solution(iPoint, 0) -1)/(c_e2-1), 2);

  F_thetat(iPoint) = min( max(FirstMaxTerm, SecondMaxTerm), 1.0);

}

void CTransLMVariable::SetF_thetat_2(unsigned long iPoint, su2double val_dist) {

  F_thetat2(iPoint) = min(F_wake(iPoint)*exp(-pow(val_dist/delta_param(iPoint), 4)), 1.0);

}

void CTransLMVariable::SetReThetat_SCF(unsigned long iPoint, su2double val_dist, su2double val_density, su2double* val_velocity, su2double val_velocityMag, su2double hRoughness, su2double* val_vorticity, su2double val_viscosity, su2double val_eddy_viscosity) {

  su2double VelocityNormalized[3];
  for(auto iDim = 0u; iDim < nDim; iDim++){
    VelocityNormalized[iDim] = val_velocity[iDim]/val_velocityMag;
  }

  su2double StreamwiseVort = 0.0;
  for(auto iDim = 0u; iDim < nDim; iDim++){
    StreamwiseVort += VelocityNormalized[iDim] * val_vorticity[iDim];
  }
  StreamwiseVort = abs(StreamwiseVort);

  su2double H_CF = StreamwiseVort * val_dist / val_velocityMag;
  su2double DeltaH_CF = H_CF * (1.0 + min(val_eddy_viscosity/val_viscosity, 0.4));
  su2double DeltaH_CF_Minus = max(-1.0*(0.1066-DeltaH_CF), 0.0);
  su2double DeltaH_CF_Plus = max(0.1066-DeltaH_CF, 0.0);
  su2double fDeltaH_CF_Minus = 75.0 * tanh(DeltaH_CF_Minus/0.0125);
  su2double fDeltaH_CF_Plus = 6200 * DeltaH_CF_Plus + 50000 * DeltaH_CF_Plus * DeltaH_CF_Plus;

  su2double toll = 1e-5;
  su2double error = toll+1.0;
  su2double rethetat_SCF_old = 20.0;   
  int nMax = 100;

  int iter;
  for (iter=0;iter<nMax && error > toll;iter++) {

    thetat_SCF(iPoint) = rethetat_SCF_old*val_viscosity/(val_density*(val_velocityMag/0.82));
    thetat_SCF(iPoint) = max(1e-20, thetat_SCF(iPoint));

    ReThetat_SCF(iPoint) = -35.088 * log(hRoughness/thetat_SCF(iPoint)) + 319.51 + fDeltaH_CF_Plus - fDeltaH_CF_Minus;

    error = abs(ReThetat_SCF(iPoint) - rethetat_SCF_old)/rethetat_SCF_old;
   
    rethetat_SCF_old = ReThetat_SCF(iPoint);

  }

}

void CTransLMVariable::Setgamma_sep(unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_viscosity,
                                    su2double val_dist, su2double StrainMag, su2double val_vort, su2double* constants) {


  FReattach(iPoint) = exp(-pow(R_T(iPoint) / 20.0, 4));

  su2double maxInnerTerm = max(0.0, (reV(iPoint) / (3.235 * rethetac(iPoint))) - 1);
  gamma_sep(iPoint) = min(constants[5] * maxInnerTerm * FReattach(iPoint), 2.0);
  gamma_sep(iPoint) = gamma_sep(iPoint) * F_thetat(iPoint);


}


