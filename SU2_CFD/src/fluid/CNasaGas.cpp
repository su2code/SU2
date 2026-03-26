/*!
 * \file CNasaGas.cpp
 * \brief Source of the NASA polynomial gas model.
 * \author Om Giri
 * \version 8.4.0 "Harrier"
 */

#include "../../include/fluid/CNasaGas.hpp"

CNasaGas::CNasaGas(su2double gamma, su2double R, bool CompEntropy) : CIdealGas(gamma, R, CompEntropy) {
  // Constructor inherited from CIdealGas
}

void CNasaGas::SetCpModel(const CConfig* config, su2double val_Temperature_Ref) {
  for (int i = 0; i < 7; ++i) {
    CpLowPolyCoefficientsND[i] = config->GetNASA_CpLowPolyCoeff(i);
    CpHighPolyCoefficientsND[i] = config->GetNASA_CpHighPolyCoeff(i);
  }
  TransitionTemperatureND = config->GetNASA_TransitionTemperature() / val_Temperature_Ref;
}

void CNasaGas::ComputeProperties_T(su2double T) {
  Temperature = T;
  const array<su2double, 7>* c;

  if (T < TransitionTemperatureND) {
    c = &CpLowPolyCoefficientsND;
  } else {
    c = &CpHighPolyCoefficientsND;
  }

  /*--- NASA polynomials:
     Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
     H/RT = -a1*T^-2 + a2*log(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + a8/T
     S/R  = -a1*T^-2/2 - a2*T^-1 + a3*log(T) + a4*T + a5*T^2/2 + a6*T^3/3 + a7*T^4/4 + a9

     Note: SU2 usually uses a 7-coefficient form:
     Cp/R = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4
     H/RT = a0 + a1*T/2 + a2*T^2/3 + a3*T^3/4 + a4*T^4/5 + a5/T
     S/R  = a0*log(T) + a1*T + a2*T^2/2 + a3*T^3/3 + a4*T^4/4 + a6

     Wait, standard NASA 7-coefficient format is:
     Cp/R = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4
     H/RT = a0 + a1*T/2 + a2*T^2/3 + a3*T^3/4 + a4*T^4/5 + a5/T
     S/R  = a0*lnT + a1*T + a2*T^2/2 + a3*T^3/3 + a4*T^4/4 + a6

     I will use this standard 7-coefficient format.
  ---*/

  Cp = ((*c)[0] + (*c)[1]*T + (*c)[2]*T*T + (*c)[3]*T*T*T + (*c)[4]*T*T*T*T) * Gas_Constant;
  Enthalpy = ((*c)[0] + (*c)[1]*T/2.0 + (*c)[2]*T*T/3.0 + (*c)[3]*T*T*T/4.0 + (*c)[4]*T*T*T*T/5.0 + (*c)[5]/T) * T * Gas_Constant;

  if (ComputeEntropy) {
    Entropy = ((*c)[0]*log(T) + (*c)[1]*T + (*c)[2]*T*T/2.0 + (*c)[3]*T*T*T/3.0 + (*c)[4]*T*T*T*T/4.0 + (*c)[6]) * Gas_Constant;
    // Note: This is partial entropy (temperature dependent part). Pressure part added in SetTDState_rhoe if needed.
    // However, SU2 CIdealGas includes log(1/rho) term.
  }

  // Cv = Cp - R
  Cv = Cp - Gas_Constant;
  // Gamma = Cp / Cv
  Gamma = Cp / Cv;
  Gamma_Minus_One = Gamma - 1.0;
}

void CNasaGas::SetTDState_rhoe(su2double rho, su2double e) {
  Density = rho;
  StaticEnergy = e;

  /*--- Solve for T using Newton-Raphson: e(T) = h(T) - R*T ---*/
  su2double T_iter = Temperature;
  if (T_iter <= 0.0) T_iter = 1.0; // Initial guess

  for (int i = 0; i < 10; ++i) {
    ComputeProperties_T(T_iter);
    su2double e_iter = Enthalpy - Gas_Constant * T_iter;
    su2double de_dT = Cv; // de/dT = Cv for ideal gas even with variable Cp
    su2double deltaT = (e - e_iter) / de_dT;
    T_iter += deltaT;
    if (abs(deltaT) < 1e-8 * T_iter) break;
  }

  Temperature = T_iter;
  ComputeProperties_T(Temperature);

  Pressure = Density * Gas_Constant * Temperature;
  SoundSpeed2 = Gamma * Pressure / Density;

  /*--- Derivatives ---*/
  dPde_rho = Density * Gas_Constant / Cv;
  dPdrho_e = Gas_Constant * Temperature - StaticEnergy * dPde_rho; // From P = rho*R*T and e = e(T)

  dTde_rho = 1.0 / Cv;
  dTdrho_e = 0.0;

  if (ComputeEntropy) {
    // Entropy was S(T). Now add pressure/density part: S = S(T) - R*ln(rho*R) + R*ln(R_ref?)
    // CIdealGas uses: Entropy = (1.0 / Gamma_Minus_One * log(Temperature) + log(1.0 / Density)) * Gas_Constant;
    // For NASA, S(T) already has the log(T) part.
    Entropy += log(1.0 / Density) * Gas_Constant;
  }
}

void CNasaGas::SetTDState_PT(su2double P, su2double T) {
  Pressure = P;
  Temperature = T;
  ComputeProperties_T(T);
  Density = P / (Gas_Constant * T);
  StaticEnergy = Enthalpy - Gas_Constant * T;
  SetTDState_rhoe(Density, StaticEnergy);
}

void CNasaGas::SetTDState_Prho(su2double P, su2double rho) {
  Pressure = P;
  Density = rho;
  Temperature = P / (Density * Gas_Constant);
  ComputeProperties_T(Temperature);
  StaticEnergy = Enthalpy - Gas_Constant * Temperature;
  SetTDState_rhoe(Density, StaticEnergy);
}

void CNasaGas::SetTDState_rhoT(su2double rho, su2double T) {
  Density = rho;
  Temperature = T;
  ComputeProperties_T(T);
  Pressure = Density * Gas_Constant * T;
  StaticEnergy = Enthalpy - Gas_Constant * T;
  SetTDState_rhoe(Density, StaticEnergy);
}

void CNasaGas::SetTDState_Ps(su2double P, su2double s) {
  /*--- This would require another nested iteration to find T from P and s.
     For now, we can use an approximate T or a simple Newton solver. ---*/
  su2double T_iter = Temperature;
  if (T_iter <= 0.0) T_iter = 1.0;

  for (int i = 0; i < 10; ++i) {
    ComputeProperties_T(T_iter);
    su2double s_iter = Entropy + log(Gas_Constant * T_iter / P) * Gas_Constant;
    su2double ds_dT = Cp / T_iter;
    su2double deltaT = (s - s_iter) / ds_dT;
    T_iter += deltaT;
    if (abs(deltaT) < 1e-8 * T_iter) break;
  }
  SetTDState_PT(P, T_iter);
}
