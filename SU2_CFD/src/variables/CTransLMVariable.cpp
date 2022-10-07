//
// Created by marcocera on 25/02/22.
//

#include "../../include/variables/CTransLMVariable.hpp"
#include <cmath>

// bisogna capire se va fatto come "if(LM)" in SST o se va fatto a parte (come questo file nuovo)

CTransLMVariable::CTransLMVariable(su2double intermittency, su2double Re_theta, const su2double* constants, unsigned long npoint, unsigned long ndim,
                                         unsigned long nvar, CConfig *config) :
                                         CTurbVariable(npoint, ndim, nvar, config) {

    for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {
        Solution(iPoint,0) = intermittency;
        Solution(iPoint,1) = Re_theta;
    }

    Solution_Old = Solution;   // queste assegnazioni a Solution e Solution_Old (Runge-Kutta problem) sono copiate da SST =>
                               // bisogna capire se cambiare nome o se il vector è pensato per contenere un numero qualsiasi
                               // di elementi (forse contiene nvar elementi)


    F_length.resize(nPoint) = su2double(0.0);
    F_onset.resize(nPoint) = su2double(0.0);
    F_turb.resize(nPoint) = su2double(0.0);
    F_onset1.resize(nPoint) = su2double(0.0);
    reV.resize(nPoint) = su2double(0.0);
    F_onset2.resize(nPoint) = su2double(0.0);
    R_T.resize(nPoint) = su2double(0.0);
    F_onset3.resize(nPoint) = su2double(0.0);
    F_length1.resize(nPoint) = su2double(0.0);
    F_sublayer.resize(nPoint) = su2double(0.0);
    rew.resize(nPoint) = su2double(0.0);
    rethetac.resize(nPoint) = su2double(0.0);
    T.resize(nPoint) = su2double(0.0);
    rethetat_eq.resize(nPoint) = su2double(0.0);
    F_thetat.resize(nPoint) = su2double(0.0);
    Velocity_Mag.resize(nPoint) = su2double(0.0);
    delta_param.resize(nPoint) = su2double(0.0);
    F_wake.resize(nPoint) = su2double(0.0);
    lambda_theta.resize(nPoint) = su2double(0.0);
    Turb_Intens.resize(nPoint) = su2double(0.0);
    du_dx.resize(nPoint) = su2double(0.0);
    du_dy.resize(nPoint) = su2double(0.0);
    du_dz.resize(nPoint) = su2double(0.0);
    dU_ds.resize(nPoint) = su2double(0.0);
    F_lambda.resize(nPoint) = su2double(0.0);
    thetat.resize(nPoint) = su2double(0.0);
    gamma_sep.resize(nPoint) = su2double(0.0);
    gamma_eff.resize(nPoint) = su2double(0.0);
    FReattach.resize(nPoint) = su2double(0.0);

    F_thetat2.resize(nPoint) = su2double(0.0);
    ReThetat_SCF.resize(nPoint) = su2double(0.0);
    thetat_SCF.resize(nPoint) = su2double(0.0);

    ProductionGamma.resize(nPoint) = su2double(0.0);
    ProductionReTheta.resize(nPoint) = su2double(0.0);
    DestructionGamma.resize(nPoint) = su2double(0.0);
    dist.resize(nPoint) = su2double(0.0);
    Ux.resize(nPoint) = su2double(0.0);
    Uy.resize(nPoint) = su2double(0.0);
    Uz.resize(nPoint) = su2double(0.0);

    LambdaTheta_L.resize(nPoint) = su2double(0.0);
    TU_L.resize(nPoint) = su2double(0.0);
    FPG.resize(nPoint) = su2double(0.0);
    ReThetaC_Corr.resize(nPoint) = su2double(0.0);
    ReThetaC_Corr_New.resize(nPoint) = su2double(0.0);
    F_Lambda_New.resize(nPoint) = su2double(0.0);

    k.resize(nPoint) = su2double(0.0);
    w.resize(nPoint) = su2double(0.0);


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


void CTransLMVariable::Set_kAndw(unsigned long iPoint, su2double VorticityMag, su2double VelocityMag, su2double val_eddy_viscosity, su2double StrainMag, su2double val_dist,
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


//    Turb_Intens(iPoint) = 100*sqrt(2*k(iPoint)/3)/VelocityMag;
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
//      cout << "Using MENTER_LANGTRY" << endl;
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
//      cout << "rethetac(" << iPoint << ") = " << rethetac(iPoint) << " ";
      break;
    }

  }

  rethetac(iPoint) = max(rethetac(iPoint), 1e-20);

}

void CTransLMVariable::FLength_Correlations(unsigned long iPoint){

  su2double ReTheta_Tilde = Solution(iPoint, 1);

  switch (Trans_Correlation) {
    case TURB_TRANS_CORRELATION::MALAN: {
//      cout << "Using MALAN" << endl;
      F_length1(iPoint) = min(exp(7.168 - 0.01173 * ReTheta_Tilde) + 0.5, 300.0);
      break;
    }

    case TURB_TRANS_CORRELATION::SULUKSNA: {
//      cout << "Using SULUKSNA" << endl;
      su2double FirstTerm = -pow(0.025 * ReTheta_Tilde, 2) + 1.47 * ReTheta_Tilde - 120.0;
      F_length1(iPoint) = min(max(FirstTerm, 125.0), ReTheta_Tilde);
      break;
    }

    case TURB_TRANS_CORRELATION::KRAUSE: {
//      cout << "Using KRAUSE" << endl;
      F_length1(iPoint) = 3.39 * ReTheta_Tilde + 55.03;
      break;
    }

    case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
//      cout << "Using MEDIDA_BAEDER" << endl;
      su2double FirstTerm = 0.171 * pow(Turb_Intens(iPoint), 2);
      su2double SecondTerm = 0.0083 * pow(Turb_Intens(iPoint), 1);
      F_length1(iPoint) = (FirstTerm - SecondTerm + 0.0306);
      break;
    }

    case TURB_TRANS_CORRELATION::MEDIDA: {
//      cout << "Using MEDIDA" << endl;
      F_length1(iPoint) = 40;
      break;
    }

    case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {
//      cout << "Using MENTER_LANGTRY" << endl;
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


//suddivisione SetQuantities in 7 funzioni minori - inizio 28/04/2022 - Marco:
void CTransLMVariable::SetF_onset(unsigned long iPoint, su2double val_density, su2double StrainMag, su2double val_dist,
                                  su2double val_viscosity, su2double *TurbVars, su2double val_eddy_viscosity, su2double VorticityMag, CConfig *config) {


  reV(iPoint) = val_density*StrainMag*val_dist*val_dist/val_viscosity;
//  reV(iPoint) = val_density*VorticityMag*val_dist*val_dist/val_viscosity;

  // If model is coupled with SST then TI can be locally computed
  if (TurbKW){
    Turb_Intens(iPoint) = 100*sqrt(2*TurbVars[0]/3)/Velocity_Mag(iPoint);
    Turb_Intens(iPoint) = max(Turb_Intens(iPoint), 0.027);
  }



  // Computing ReTheta_C from correlations
  // For LM-SST couple, the preferred one is LANGTRY_MENTER
  // For LM-SA couple, the preferred one is MALAN
  // Change these in the configuration file!
  ReThetaC_Correlations(iPoint);



  //    cout << "Dopo rethetac" << endl;

  /*--- F_onset ---*/
  F_onset1(iPoint) = reV(iPoint)/(2.193*rethetac(iPoint));

  // If model is coupled with SA then R_T is expressed differently
  if (TurbSA)
    R_T(iPoint) = val_eddy_viscosity / val_viscosity;

  // Original formulation
  if (TurbKW)
    R_T(iPoint) = val_density*TurbVars[0]/(val_viscosity*TurbVars[1]);

  // Cambiata da me
//  if (TurbKW)
//    R_T(iPoint) = val_eddy_viscosity / val_viscosity;

  if (TurbSA && config->GetConvertSA2SST())
    R_T(iPoint) = val_density*k(iPoint)/(val_viscosity*w(iPoint));

//  cout << R_T(iPoint) << endl;


  if (TurbSA) {
    F_onset2(iPoint) = min(max(F_onset1(iPoint), pow(F_onset1(iPoint), 4)), 4.0);
    F_onset3(iPoint) = max(2.0 - pow(0.4 * R_T(iPoint), 3), 0.0);

//    F_onset2(iPoint) = min(max(F_onset1(iPoint), pow(F_onset1(iPoint), 4)), 2.0);
//    F_onset3(iPoint) = max(1.0 - pow(0.4 * R_T(iPoint), 3), 0.0);
  }
  if (TurbKW) {
    F_onset2(iPoint) = min(max(F_onset1(iPoint), pow(F_onset1(iPoint), 4)), 2.0);
    F_onset3(iPoint) = max(1.0 - pow(0.4 * R_T(iPoint), 3), 0.0);
  }

  F_onset(iPoint) = max(F_onset2(iPoint)-F_onset3(iPoint),0.0);

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
    //    cout << "Dopo F_sublayer" << endl;
    F_length(iPoint) = F_length1(iPoint) * (1 - F_sublayer(iPoint)) + 40.0 * F_sublayer(iPoint);

  }

  if(TurbSA && config->GetConvertSA2SST()) {

//    rew(iPoint) = val_density * w(iPoint) * val_dist * val_dist / val_viscosity;  // assunto "val_viscosity" come la laminar
//
//
//    F_sublayer(iPoint) = exp(-pow(0.005 * rew(iPoint), 2));
//    //    cout << "Dopo F_sublayer" << endl;
//    F_length(iPoint) = F_length1(iPoint) * (1 - F_sublayer(iPoint)) + 40.0 * F_sublayer(iPoint);

    F_length(iPoint) = F_length1(iPoint);

  }

}

void CTransLMVariable::SetF_turb(unsigned long iPoint) {

  F_turb(iPoint) = exp(-pow(0.25*R_T(iPoint),4));

}

void CTransLMVariable::Setrethetat_eq(unsigned long iPoint, su2double *Velocity, su2double VelocityMag, CMatrixView<const su2double> Velocity_Gradient,
                                      su2double *TurbVars, su2double val_viscosity, su2double val_density, CConfig *config) {


  Ux(iPoint) = Velocity[0];
  Uy(iPoint) = Velocity[1];
  if(nDim == 3) Uz(iPoint) = Velocity[2];

  Velocity_Mag(iPoint) = VelocityMag;

  // Aggiunto da me
  // Calcolate come nel paper
  du_dx(iPoint) = 0.0;
  for(auto iDim = 0u; iDim < nDim; iDim++){
    du_dx(iPoint) += 2 * Velocity[iDim] * Velocity_Gradient[iDim][0];
  }
  du_dx(iPoint) = du_dx(iPoint) * 0.5 / Velocity_Mag(iPoint);

  du_dy(iPoint) = 0.0;
  for(auto iDim = 0u; iDim < nDim; iDim++){
    du_dy(iPoint) += 2 * Velocity[iDim] * Velocity_Gradient[iDim][1];
  }
  du_dy(iPoint) = du_dy(iPoint) * 0.5 / Velocity_Mag(iPoint);

  if(nDim == 3) {
    du_dz(iPoint) = 0.0;
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      du_dz(iPoint) += 2 * Velocity[iDim] * Velocity_Gradient[iDim][2];
    }
    du_dz(iPoint) = du_dz(iPoint) * 0.5 / Velocity_Mag(iPoint);
  }
  //    cout << "Dopo du_dz" << endl;



  dU_ds(iPoint) = Velocity[0] * du_dx(iPoint) / Velocity_Mag(iPoint) + Velocity[1] * du_dy(iPoint) / Velocity_Mag(iPoint);
  if (nDim==3) {
    dU_ds(iPoint) += Velocity[2] * du_dz(iPoint) / Velocity_Mag(iPoint);
  }
  //    cout << "Dopo dU_ds" << endl;

  // If model is coupled with SA then TI is the freestream one



  //    cout << "Dopo F_Length" << endl;

  // Guess iniziale fatta con lambda = 0.
  F_lambda(iPoint) = 1.0;
  su2double toll = 1e-5;
  su2double error = toll+1.0;
  su2double rethetat_eq_old = 20.0;   // ho messo 20.0, dato che è il limite minimo imposto nell'articolo NASA, prima
                                     // era 1.0
  int nMax = 100;
  //    cout << "Dopo F_lambda" << endl;

  // Aggiustare assolutamente questa soluzione
  int iter;
  for (iter=0;iter<nMax && error > toll;iter++) {   // quante iterazioni? non è meglio un while che abbia una tolleranza sensata? nel caso sì, quale valore?

    thetat(iPoint) = rethetat_eq_old*val_viscosity/(val_density*Velocity_Mag(iPoint));
    lambda_theta(iPoint) = val_density*thetat(iPoint)*thetat(iPoint)*dU_ds(iPoint)/val_viscosity;
    // lambda_theta dovrebbe essere limitata tra -0.1 e 0.1
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

    rethetat_eq(iPoint) = max(20.0, rethetat_eq(iPoint));   // limite inferiore preso da pag 4 del sito NASA

    error = abs(rethetat_eq(iPoint) - rethetat_eq_old)/rethetat_eq_old;
    //        if(iPoint == 20){
    //          cout << "iter = " << iter << " Re_theta_t_eq = " << rethetat_eq(iPoint) << " error = " << error  << endl;
    //        }

    rethetat_eq_old = rethetat_eq(iPoint);


    }

    if(iter >= nMax) cout << "cazzo" << endl;

    // Computation of Re_theta_c with only  correlations
    su2double LambdaTheta_L_Here = 0.0;
    LambdaTheta_L_Here = -7.57e-3 * dU_ds(iPoint) * dist(iPoint)*dist(iPoint) / (val_viscosity/val_density) + 0.0128;
    LambdaTheta_L_Here = min( max(LambdaTheta_L_Here, -1.0), 1.0);

    LambdaTheta_L(iPoint) = LambdaTheta_L_Here;

    if (LambdaTheta_L(iPoint)<=0.0) {
      su2double FirstTerm = 12.986*LambdaTheta_L(iPoint);
      su2double SecondTerm = 123.66*pow(LambdaTheta_L(iPoint),2);
      su2double ThirdTerm = 405.689*pow(LambdaTheta_L(iPoint),3);
      F_Lambda_New(iPoint) = 1 + ( FirstTerm + SecondTerm + ThirdTerm ) * exp(- pow(Turb_Intens(iPoint)/1.5,1.5));
    }
    else {
      su2double FirstTerm = exp(-35.0*LambdaTheta_L(iPoint));
      F_Lambda_New(iPoint) = 1.0 + 0.275*(1.0 - FirstTerm) * exp(-2.0*Turb_Intens(iPoint));
    }

    if (Turb_Intens(iPoint)<=1.3) {
      su2double FirstTerm = 589.428*Turb_Intens(iPoint);
      su2double SecondTerm = 0.2196/(Turb_Intens(iPoint)*Turb_Intens(iPoint));
      ReThetaC_Corr_New(iPoint) = (1173.51 - FirstTerm + SecondTerm) * F_Lambda_New(iPoint);
    }
    else {
      ReThetaC_Corr_New(iPoint) = 331.5*pow(Turb_Intens(iPoint)-0.5658,-0.671)*F_Lambda_New(iPoint);
    }

    ReThetaC_Corr_New(iPoint) = max(20.0, ReThetaC_Corr_New(iPoint));   // limite inferiore preso da pag 4 del sito NASA


    if (TurbKW){
      // Tentativo per one equation transition model
      // Taken from "A One-Equation Local Correlation-Based Transition Model" DOI 10.1007/s10494-015-9622-4
      su2double TU_L_Here = 0.0;
      TU_L_Here = min(100*sqrt(2*TurbVars[0]/3)/(TurbVars[1]*dist(iPoint)), 100.0 );
      TU_L(iPoint) = TU_L_Here;

      su2double CPG1 = 14.68;
      su2double CPG2 = -7.34;
      su2double CPG3 = 0.0;
      su2double CPG1_lim = 1.5;
      su2double CPG2_lim = 3.0;

      su2double FPG_Here = 0.0;
      if(LambdaTheta_L_Here >= 0){

        FPG_Here = min(1 + CPG1*LambdaTheta_L_Here, CPG1_lim);

      }
      else {

        FPG_Here = min(1 + CPG2*LambdaTheta_L_Here + CPG3*min(LambdaTheta_L_Here+0.0681, 0.0), CPG2_lim);

      }

      FPG_Here = max(FPG_Here, 0.0);
      FPG(iPoint) = FPG_Here;

      su2double CTU1 = 100.0;
      su2double CTU2 = 1000.0;
      su2double CTU3 = 1.0;
      su2double ReThetaC_Here = 0.0;

      ReThetaC_Here = CTU1 + CTU2*exp(-CTU3*TU_L_Here*FPG_Here);
      ReThetaC_Corr(iPoint) = ReThetaC_Here;
    }





  //    cout << "iter = " << iter << " Re_theta_t_eq(" << iPoint << ") = " << rethetat_eq(iPoint) << " error = " << error  << endl;

  //    cout << "iter = " << iter << endl;

}

void CTransLMVariable::SetT(unsigned long iPoint, su2double *Velocity, su2double val_viscosity, su2double val_density, su2double val_localGridLength, su2double val_mu_T) {


  T(iPoint) = 500*val_viscosity/(val_density*Velocity_Mag(iPoint)*Velocity_Mag(iPoint));   // assunto "val_viscosity" come la laminar

  if(TransModel == TURB_TRANS_MODEL::LM2015){
    T(iPoint) = min(T(iPoint), val_density*val_localGridLength*val_localGridLength/(val_viscosity + val_mu_T));
  }

  T(iPoint) = max(T(iPoint), 1e-20);

}

void CTransLMVariable::SetF_thetat(unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_dist,
                                   su2double val_viscosity, su2double val_vort, su2double* constants, CConfig *config) {


  su2double theta_BL = Solution(iPoint, 1) * val_viscosity / (val_density * Velocity_Mag(iPoint));
  su2double delta_BL = 15.0 * theta_BL / 2.0;
  delta_param(iPoint) = 50.0 * val_vort * val_dist * delta_BL / Velocity_Mag(iPoint);
  delta_param(iPoint) = max(delta_param(iPoint), 1e-20);


  // Forse questa distinzione non serve visto che rew viene inizializzato a 0
  // If model is coupled with SA then F_wake = 1
  if (TurbSA)
    F_wake(iPoint) = 1.0;

  if (TurbKW)
    F_wake(iPoint) = exp(-pow(0.00001*rew(iPoint),2));

  if (TurbSA && config->GetConvertSA2SST())
    F_wake(iPoint) = 1.0;
//    F_wake(iPoint) = exp(-pow(0.00001*rew(iPoint),2));


  su2double c_e2 = constants[3];
  su2double FirstMaxTerm = F_wake(iPoint) * exp(-pow(val_dist/delta_param(iPoint), 4));
  su2double SecondMaxTerm = 1 - pow((c_e2 * Solution(iPoint, 0) -1)/(c_e2-1), 2);
  //    cout << "Dopo SecondMaxTerm" << endl;

  F_thetat(iPoint) = min( max(FirstMaxTerm, SecondMaxTerm), 1.0);
  //    cout << "F_thetat(" << iPoint << ") = " << F_thetat(iPoint) << endl;

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

  // Guess iniziale fatta con lambda = 0.
  su2double toll = 1e-5;
  su2double error = toll+1.0;
  su2double rethetat_SCF_old = 20.0;   // ho messo 20.0, dato che non sapevo che mettere
  int nMax = 100;
  //    cout << "Dopo F_lambda" << endl;

  int iter;
  for (iter=0;iter<nMax && error > toll;iter++) {

    thetat_SCF(iPoint) = rethetat_SCF_old*val_viscosity/(val_density*(val_velocityMag/0.82));
    thetat_SCF(iPoint) = max(1e-20, thetat_SCF(iPoint));

    ReThetat_SCF(iPoint) = -35.088 * log(hRoughness/thetat_SCF(iPoint)) + 319.51 + fDeltaH_CF_Plus - fDeltaH_CF_Minus;

//    ReThetat_SCF(iPoint) = max(20.0, ReThetat_SCF(iPoint));   // limite inferiore preso da pag 4 del sito NASA

    error = abs(ReThetat_SCF(iPoint) - rethetat_SCF_old)/rethetat_SCF_old;
    //        if(iPoint == 20){
    //          cout << "iter = " << iter << " Re_theta_t_eq = " << rethetat_eq(iPoint) << " error = " << error  << endl;
    //        }

    rethetat_SCF_old = ReThetat_SCF(iPoint);


  }

  if(iter >= nMax) cout << "cazzo" << endl;

}

void CTransLMVariable::Setgamma_sep(unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_viscosity,
                                    su2double val_dist, su2double StrainMag, su2double val_vort, su2double* constants) {


  FReattach(iPoint) = exp(-pow(R_T(iPoint) / 20.0, 4));

  su2double maxInnerTerm = max(0.0, (reV(iPoint) / (3.235 * rethetac(iPoint))) - 1);
  gamma_sep(iPoint) = min(constants[5] * maxInnerTerm * FReattach(iPoint), 2.0);
  gamma_sep(iPoint) = gamma_sep(iPoint) * F_thetat(iPoint);


  //    cout << " gamma_sep(" << iPoint << ") = " << gamma_sep(iPoint) << endl;
}



void CTransLMVariable::SetDebugQuantities(unsigned long iPoint, su2double* constants, su2double val_density,
                                          su2double StrainMag, su2double val_vort, su2double val_dist,
                                          su2double val_viscosity, su2double *TurbVars) {

  ProductionGamma(iPoint) = F_length(iPoint) * constants[0] * val_density * StrainMag * pow(Solution(iPoint, 0) * F_onset(iPoint), 0.5) *
                            (1.0 - constants[2] * Solution(iPoint, 0));

  ProductionReTheta(iPoint) = constants[4] * (val_density / T(iPoint)) * (rethetat_eq(iPoint) - Solution(iPoint, 1)) * (1.0 - F_thetat(iPoint));
  DestructionGamma(iPoint) = constants[1] * val_density * val_vort * Solution(iPoint, 0) * F_turb(iPoint) * (constants[3] * Solution(iPoint, 0) - 1.0);

  dist(iPoint) = val_dist;

  //    cout << "ProductionReTheta(" << iPoint << ") = " <<ProductionReTheta(iPoint) << endl;

  bool print = false;
  if(print) {
    int iPoint2print = 200;
    if (iPoint == iPoint2print) {
      cout << "F_length(0) = " << F_length(iPoint2print) << endl;
      cout << "F_onset(0) = " << F_onset(iPoint2print) << endl;
      cout << "F_turb(0) = " << F_turb(iPoint2print) << endl;
      cout << "val_density = " << val_density << endl;
      cout << "val_viscosity = " << val_viscosity << endl;
      cout << "TurbVars(0) = " << TurbVars[0] << endl;
      cout << "TurbVars(1) = " << TurbVars[1] << endl;
      cout << "F_onset1(0) = " << F_onset1(iPoint2print) << endl;
      cout << "reV(0) = " << reV(iPoint2print) << endl;
      cout << "F_onset2(0) = " << F_onset2(iPoint2print) << endl;
      cout << "R_T(0) = " << R_T(iPoint2print) << endl;
      cout << "F_onset3(0) = " << F_onset3(iPoint2print) << endl;
      cout << "F_length1(0) = " << F_length1(iPoint2print) << endl;
      cout << "F_sublayer(0) = " << F_sublayer(iPoint2print) << endl;
      cout << "rew(0) = " << rew(iPoint2print) << endl;
      cout << "rethetac(0) = " << rethetac(iPoint2print) << endl;
      cout << "T(0) = " << T(iPoint2print) << endl;
      cout << "rethetat_eq(0) = " << rethetat_eq(iPoint2print) << endl;
      cout << "F_thetat(0) = " << F_thetat(iPoint2print) << endl;
      cout << "Velocity_Mag(0) = " << Velocity_Mag(iPoint2print) << endl;
      cout << "delta_param(0) = " << delta_param(iPoint2print) << endl;
      cout << "F_wake(0) = " << F_wake(0) << endl;
      cout << "lambda_theta(0) = " << lambda_theta(iPoint2print) << endl;
      cout << "Turb_Intens(0) = " << Turb_Intens(iPoint2print) << endl;
      cout << "du_dx(0) = " << du_dx(iPoint2print) << endl;
      cout << "du_dy(0) = " << du_dy(iPoint2print) << endl;
      cout << "du_dz(0) = " << du_dz(iPoint2print) << endl;
      cout << "dU_ds(0) = " << dU_ds(iPoint2print) << endl;
      cout << "F_lambda(0) = " << F_lambda(iPoint2print) << endl;
      cout << "thetat(0) = " << thetat(iPoint2print) << endl;
      cout << "gamma_sep(0) = " << gamma_sep(iPoint2print) << endl;
      cout << "gamma_eff(0) = " << gamma_eff(iPoint2print) << endl;
    }
  }
}
