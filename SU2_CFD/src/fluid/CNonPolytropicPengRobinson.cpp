/*!
 * \file CNonPolytropicPengRobinson.cpp
 * \brief Source of the non-polytropic Peng-Robinson model.
 * \author B. Fuentes Monjas
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

#include "../../include/fluid/CNonPolytropicPengRobinson.hpp"

CNonPolytropicPengRobinson::CNonPolytropicPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w, const CConfig *config)
    : CPengRobinson (gamma, R, Pstar, Tstar, w) {
  Temperature = TstarCrit; // Initialization of Temperature to critical value
}

void CNonPolytropicPengRobinson::SetCpModel(const CConfig* config) {
  for (int i = 0; i < config->GetnPolyCoeffs(); ++i) {
    coeffs_[i] = config->GetCp_PolyCoeffND(i);
  }
}

su2double CNonPolytropicPengRobinson::T_v_h(su2double v, su2double h){
  su2double toll = 1e-9, FACTOR = 0.2;
  su2double f, fv, C0, C1, C2, C4, C6, C8, C10, T, d, atanh;
  su2double x1, x2, xmid, dx, fx1, fx2, fmid = 1.0, rtb;
  su2double sqrt2 = sqrt(2.0);
  unsigned short countrtb = 0, NTRY = 100, ITMAX = 100;

  d = (v * v + 2 * b * v - b * b);

  atanh = (log(1.0 + (b * sqrt2 / (v + b))) - log(1.0 - (b * sqrt2 / (v + b)))) / 2.0;

  fv = atanh;

  C0 = -(h + a * (1 + k) * (1 + k) * (fv / (b * sqrt2) + v / d));
  C1 = a * k * (k + 1) / sqrt(TstarCrit) * (fv / (b * sqrt2) + 2 * v / d);
  C2 = coeffs_[0] - Gas_Constant + (Gas_Constant * v / (v - b)) - a * v * k * k / (TstarCrit * d);
  C4 = coeffs_[1] / 2;
  C6 = coeffs_[2] / 3;
  C8 = coeffs_[3] / 4;
  C10 = coeffs_[4] / 5;

  T = h / ComputeIntegralCp0_DT(Temperature) * Temperature;

  x1 = sqrt(0.9*T);
  x2 = sqrt(1.1*T);

  fx1 = C0 + (C1 * x1) + (C2 * x1 * x1) + (C4 * x1 * x1 * x1 * x1) + (C6 * x1 * x1 * x1 * x1 * x1 * x1) +
        (C8 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1) +
        (C10 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1);
  fx2 = C0 + (C1 * x2) + (C2 * x2 * x2) + (C4 * x2 * x2 * x2 * x2) + (C6 * x2 * x2 * x2 * x2 * x2 * x2) +
        (C8 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2) + 
        (C10 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2);

  // zbrac algorithm NR

  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        fx1 = C0 + (C1 * x1) + (C2 * x1 * x1) + (C4 * x1 * x1 * x1 * x1) + (C6 * x1 * x1 * x1 * x1 * x1 * x1) +
              (C8 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1) + 
              (C10 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1);
      } else {
        x2 += FACTOR * (x2 - x1);
        fx2 = C0 + (C1 * x2) + (C2 * x2 * x2) + (C4 * x2 * x2 * x2 * x2) + (C6 * x2 * x2 * x2 * x2 * x2 * x2) +
              (C8 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2) +
              (C10 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2);
      }
    }
  }

  // rtbis algorithm NR

  f = fx1;
  fmid = fx2;
  if (f * fmid >= 0.0) {
    cout << "Root must be bracketed for bisection in rtbis\n"
            "No solution in the temperature range specified" << endl;
  }
  rtb = f < 0.0 ? (static_cast<void>(dx = x2 - x1), x1) : (static_cast<void>(dx = x1 - x2), x2);
  do {
    xmid = rtb + (dx *= 0.5);
    fmid = C0 + (C1 * xmid) + (C2 * xmid * xmid) + (C4 * xmid * xmid * xmid * xmid) + 
           (C6 * xmid * xmid * xmid * xmid * xmid * xmid) +
           (C8 * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid) +
           (C10 * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid);

    if (fmid <= 0.0) rtb = xmid;
    countrtb++;
  } while (abs(fmid) > toll && countrtb < ITMAX);

  T = xmid;

  return T * T;
}

void CNonPolytropicPengRobinson::SetTDState_rhoe(su2double rho, su2double e){
  su2double toll = 1e-9, FACTOR = 0.2;
  su2double f, fv, A, B, C0, C1, C2, C4, C6, C8, C10, T, atanh;
  su2double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
  unsigned short countrtb = 0, NTRY = 100, ITMAX = 100;
  su2double DpDd_T, DpDT_d, DeDd_T, Cv0;
  su2double sqrt2, a2T, rho2, sig, sres;

  Density = rho;
  StaticEnergy = e;

  AD::StartPreacc();
  AD::SetPreaccIn(rho);
  AD::SetPreaccIn(e);

  rho2 = rho * rho;
  sqrt2 = sqrt(2.0);

  atanh = (log(1.0 + (rho * b * sqrt2 / (1 + rho * b))) - log(1.0 - (rho * b * sqrt2 / (1 + rho * b)))) / 2.0;

  fv = atanh;

  C0 = -(a * (k + 1) * (k + 1) * fv / (b * sqrt2) + e);
  C1 = a * k * (k + 1) * fv / (b * sqrt2 * sqrt(TstarCrit));
  C2 = coeffs_[0] - Gas_Constant;
  C4 = coeffs_[1] / 2;
  C6 = coeffs_[2] / 3;
  C8 = coeffs_[3] / 4;
  C10 = coeffs_[4] / 5;

  T = e / (ComputeIntegralCp0_DT(Temperature) - Gas_Constant*Temperature) * Temperature;

  x1 = sqrt(0.9*T);
  x2 = sqrt(1.1*T);

  fx1 = C0 + (C1 * x1) + (C2 * x1 * x1) + (C4 * x1 * x1 * x1 * x1) + (C6 * x1 * x1 * x1 * x1 * x1 * x1) +
        (C8 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1) +
        (C10 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1);
  fx2 = C0 + (C1 * x2) + (C2 * x2 * x2) + (C4 * x2 * x2 * x2 * x2) + (C6 * x2 * x2 * x2 * x2 * x2 * x2) +
        (C8 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2) +
        (C10 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2);
  
  // cout << "Intervals: " << fx1 << " and " << fx2 << endl;

  // zbrac algorithm NR

  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        fx1 = C0 + (C1 * x1) + (C2 * x1 * x1) + (C4 * x1 * x1 * x1 * x1) + (C6 * x1 * x1 * x1 * x1 * x1 * x1) +
              (C8 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1) +
              (C10 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1 * x1);
      } else {
        x2 += FACTOR * (x2 - x1);
        fx2 = C0 + (C1 * x2) + (C2 * x2 * x2) + (C4 * x2 * x2 * x2 * x2) + (C6 * x2 * x2 * x2 * x2 * x2 * x2) +
              (C8 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2) +
              (C10 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2);
      }
    }
  }

  // rtbis algorithm NR

  // cout << "Intervals2: " << fx1 << " and " << fx2 << endl;
  f = fx1;
  fmid = fx2;
  if (f * fmid >= 0.0) {
    cout << "Root must be bracketed for bisection in rtbis in the main thermodynamics routine" << endl;
  }
  rtb = f < 0.0 ? (static_cast<void>(dx = x2 - x1), x1) : (static_cast<void>(dx = x1 - x2), x2);
  do {
    xmid = rtb + (dx *= 0.5);
    fmid = C0 + (C1 * xmid) + (C2 * xmid * xmid) + (C4 * xmid * xmid * xmid * xmid) + 
           (C6 * xmid * xmid * xmid * xmid * xmid * xmid) +
           (C8 * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid) +
           (C10 * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid * xmid);

    if (fmid <= 0.0) rtb = xmid;
    countrtb++;
  } while (abs(fmid) > toll && countrtb < ITMAX);

  T = xmid;
  
  Temperature = T*T;

  a2T = alpha2(Temperature);

  A = (1 / rho2 + 2 * b / rho - b * b);
  B = 1 / rho - b;

  Pressure = Temperature * Gas_Constant / B - a * a2T / A;

  sig = ComputeIntegralCp0_T_DT(Temperature) - Gas_Constant * log(Temperature * Density);
  sres = - Gas_Constant * log(1 / (rho * B)) -
         a * sqrt(a2T) * k * fv / (b * sqrt2 * sqrt(Temperature * TstarCrit));

  Entropy = sig + sres;

  DpDd_T = (Temperature * Gas_Constant / (B * B) - 2 * a * a2T * (1 / rho + b) / (A * A)) / (rho2);

  DpDT_d = Gas_Constant / B + a * k / A * sqrt(a2T / (Temperature * TstarCrit));

  Cv0 = GetCv0();

  Cv = Cv0 + (a * k * (k + 1) * fv) / (2 * b * sqrt(2 * Temperature * TstarCrit));

  dPde_rho = DpDT_d / Cv;

  DeDd_T = -a * (1 + k) * sqrt(a2T) / A / (rho2);

  dPdrho_e = DpDd_T - dPde_rho * DeDd_T;

  SoundSpeed2 = dPdrho_e + Pressure / (rho2)*dPde_rho;

  dTde_rho = 1 / Cv;

  Zed = Pressure / (Gas_Constant * Temperature * rho);

  AD::SetPreaccOut(Temperature);
  AD::SetPreaccOut(SoundSpeed2);
  AD::SetPreaccOut(dPde_rho);
  AD::SetPreaccOut(dPdrho_e);
  AD::SetPreaccOut(Zed);
  AD::SetPreaccOut(dTde_rho);
  AD::SetPreaccOut(Pressure);
  AD::SetPreaccOut(Entropy);
  AD::EndPreacc();
}

void CNonPolytropicPengRobinson::SetTDState_PT(su2double P, su2double T){
  su2double toll = 1e-6;
  su2double A, B, Z, DZ = 1.0, F, F1, atanh;
  su2double rho, fv, e;
  su2double sqrt2 = sqrt(2.0);
  unsigned short nmax = 20, count = 0;

  AD::StartPreacc();
  AD::SetPreaccIn(P);
  AD::SetPreaccIn(T);

  A = a * alpha2(T) * P / (T * Gas_Constant) / (T * Gas_Constant);
  B = b * P / (T * Gas_Constant);

  if (Zed > 0.1)
    Z = min(Zed, 0.99);
  else
    Z = 0.99;

  do {
    F = Z * Z * Z + Z * Z * (B - 1.0) + Z * (A - 2 * B - 3 * B * B) + (B * B * B + B * B - A * B);
    F1 = 3 * Z * Z + 2 * Z * (B - 1.0) + (A - 2 * B - 3 * B * B);
    DZ = F / F1;
    Z -= DZ;
  } while (abs(DZ) > toll && count < nmax);

  if (count == nmax) {
    cout << "Warning Newton-Raphson exceed number of max iteration in PT" << endl;
    cout << "Compressibility factor  " << Z << " would be substituted with " << Zed << endl;
  }
  // check if the solution is physical otherwise uses previous point  solution
  if (Z <= 1.0001 && Z >= 0.05 && count < nmax) Zed = Z;

  rho = P / (Zed * Gas_Constant * T);

  atanh = (log(1.0 + (rho * b * sqrt2 / (1 + rho * b))) - log(1.0 - (rho * b * sqrt2 / (1 + rho * b)))) / 2.0;

  fv = atanh;

  e = ComputeIntegralCp0_DT(T) - Gas_Constant * T - a * (k + 1) * sqrt(alpha2(T)) * fv / (b * sqrt2);

  AD::SetPreaccOut(rho);
  AD::SetPreaccOut(e);
  AD::EndPreacc();

  SetTDState_rhoe(rho, e);
}

void CNonPolytropicPengRobinson::SetTDState_hs(su2double h, su2double s) {
  su2double T, fv, sqrt2 = sqrt(2.0);
  su2double f, v, atanh, sig, sres;
  su2double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
  su2double toll = 1e-9, FACTOR = 0.2;
  su2double cons_s, cons_h;
  unsigned short countrtb = 0, NTRY = 100, ITMAX = 100;

  T = h / ComputeIntegralCp0_DT(Temperature) * Temperature;
  v = exp((-ComputeIntegralCp0_T_DT(T) + s)/ Gas_Constant + log(T));

  x1 = 0.9*v;
  x2 = 1.1*v;
  
  T = T_v_h(x1, h);
  atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
  fv = atanh;

  sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x1);
  sres = - Gas_Constant*log((x1)/(x1 - b)) -
         a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

  fx1 = sig + sres - s;

  T = T_v_h(x2, h);

  atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
  fv = atanh;

  sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x2);
  sres = - Gas_Constant*log((x2)/(x2 - b)) -
         a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

  fx2 = sig + sres - s;

  // zbrac algorithm NR
  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        T = T_v_h(x1, h);
        atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
        fv = atanh;
        sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x1);
        sres = - Gas_Constant*log((x1)/(x1 - b)) -
               a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));
        fx1 = sig + sres - s;
        } else {
        x2 += FACTOR * (x2 - x1);
        T = T_v_h(x2, h);
        atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
        fv = atanh;
        sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x2);
        sres = - Gas_Constant*log((x2)/(x2 - b)) -
               a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));
        fx2 = sig + sres - s;
      }
    }
  }

  // rtbis algorithm NR
  f = fx1;
  fmid = fx2;
  if (f * fmid >= 0.0) {
    cout << "Root must be bracketed for bisection in rtbis" << endl;
    SetTDState_rhoT(Density, Temperature);
  }
  rtb = f < 0.0 ? (static_cast<void>(dx = x2 - x1), x1) : (static_cast<void>(dx = x1 - x2), x2);
  do {
    xmid = rtb + (dx *= 0.5);
    T = T_v_h(xmid, h);
    atanh = (log(1.0 + (b * sqrt2 / (xmid + b))) - log(1.0 - (b * sqrt2 / (xmid + b)))) / 2.0;
    fv = atanh;
    sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / xmid);
    sres = - Gas_Constant*log((xmid)/(xmid - b)) -
            a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));
    fmid = sig + sres - s;
    if (fmid <= 0.0) rtb = xmid;
    countrtb++;
  } while (abs(fmid) > toll && countrtb < ITMAX);

  v = xmid;
  if (countrtb == ITMAX) {
    cout << "Too many bisections in rtbis" << endl;
    cout << countrtb << endl;
  }
  if (v != v) {
    cout << "not physical solution found, h and s input " << h << " " << s << endl;
    SetTDState_rhoT(Density, Temperature);
  }

  T = T_v_h(v, h);
  SetTDState_rhoT(1 / v, T);

  // consistency check
  cons_h = abs(((StaticEnergy + Pressure / Density) - h) / h);
  cons_s = abs((Entropy - s) / s);

  if (cons_h > 1e-4 || cons_s > 1e-4) {
    cout << "TD consistency not verified in hs call" << endl;
  }
}

void CNonPolytropicPengRobinson::SetEnergy_Prho(su2double P, su2double rho) {
  su2double ad;
  su2double A, B, C, T, vb1, vb2, atanh;

  AD::StartPreacc();
  AD::SetPreaccIn(P);
  AD::SetPreaccIn(rho);

  vb1 = (1 / rho - b);
  vb2 = (1 / rho / rho + 2 * b / rho - b * b);

  A = Gas_Constant / vb1 - a * k * k / TstarCrit / vb2;

  B = 2 * a * k * (k + 1) / sqrt(TstarCrit) / vb2;

  C = -P - a * (1 + k) * (1 + k) / vb2;

  T = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
  T *= T;

  atanh = (log(1.0 + (rho * b * sqrt(2.0) / (1 + rho * b))) - log(1.0 - (rho * b * sqrt(2.0) / (1 + rho * b)))) / 2.0;
  ad = a * (k + 1) * sqrt(alpha2(T)) / (b * sqrt(2.0)) * atanh;

  StaticEnergy = ComputeIntegralCp0_DT(T) - Gas_Constant * T - ad;

  AD::SetPreaccOut(StaticEnergy);
  AD::EndPreacc();
}

void CNonPolytropicPengRobinson::SetTDState_rhoT(su2double rho, su2double T) {
  su2double fv, e, atanh;

  atanh = (log(1.0 + (rho * b * sqrt(2.0) / (1 + rho * b))) - log(1.0 - (rho * b * sqrt(2.0) / (1 + rho * b)))) / 2.0;
  fv = atanh;
  e = ComputeIntegralCp0_DT(T) - Gas_Constant * T - a * (k + 1) * sqrt(alpha2(T)) / (b * sqrt(2.0)) * fv;
  SetTDState_rhoe(rho, e);
}


void CNonPolytropicPengRobinson::SetTDState_Ps(su2double P, su2double s) {
  su2double T, rho, v, cons_P, cons_s, fv, atanh, sig, sres;
  su2double x1, x2, fx1, fx2, f, fmid, rtb, dx, xmid, sqrt2 = sqrt(2.0);
  su2double toll = 1e-5, FACTOR = 0.2;
  unsigned short count = 0, NTRY = 100, ITMAX = 100;

  T = exp((s + Gas_Constant*log(P))/(ComputeIntegralCp0_T_DT(Temperature))*log(Temperature));
  v = (T * Gas_Constant) / P;

  if (Zed < 0.9999) {
    x1 = Zed * v;
    x2 = v;

  } else {
    x1 = 0.2 * v;
    x2 = v;
  }
  T = T_P_rho(P, 1.0 / x1);

  atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
  fv = atanh;
  
  sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x1);
  sres = - Gas_Constant + log((x1)/(x1 - b)) -
         a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

  fx1 = sig + sres - s;
  
  T = T_P_rho(P, 1.0 / x2);

  atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
  fv = atanh;

  sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x2);
  sres = - Gas_Constant + log((x2)/(x2 - b)) -
         a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

  fx2 = sig + sres - s;

  // zbrac algorithm NR

  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        T = T_P_rho(P, 1.0 / x1);

        atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
        fv = atanh;
        
        sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x1);
        sres = - Gas_Constant + log((x1)/(x1- b)) -
               a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

        fx1 = sig + sres - s;      
        } else {
        T = T_P_rho(P, 1.0 / x2);

        atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
        fv = atanh;
        
        sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / x2);
        sres = - Gas_Constant + log((x2)/(x2- b)) -
               a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

        fx2 = sig + sres - s;  
      }
    }
  }

  // rtbis algorithm NR

  f = fx1;
  fmid = fx2;
  if (f * fmid >= 0.0) {
    cout << "Root must be bracketed for bisection in rtbis" << endl;
    SetTDState_rhoT(Density, Temperature);
  }
  rtb = f < 0.0 ? (static_cast<void>(dx = x2 - x1), x1) : (static_cast<void>(dx = x1 - x2), x2);
  do {
    xmid = rtb + (dx *= 0.5);
    T = T_P_rho(P, 1.0 / xmid);

    atanh = (log(1.0 + (b * sqrt2 / (xmid + b))) - log(1.0 - (b * sqrt2 / (xmid + b)))) / 2.0;
    fv = atanh;
        
    sig = ComputeIntegralCp0_T_DT(T) - Gas_Constant * log(T / xmid);
    sres = - Gas_Constant + log((xmid)/(xmid- b)) -
           a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit));

    fmid = sig + sres - s;

    if (fmid <= 0.0) rtb = xmid;
    count++;
  } while (abs(fmid) > toll && count < ITMAX);

  if (count == ITMAX) {
    cout << "Too many bisections in rtbis" << endl;
  }

  rho = 1.0 / xmid;
  T = T_P_rho(P, rho);
  SetTDState_rhoT(rho, T);

  cons_P = abs((Pressure - P) / P);
  cons_s = abs((Entropy - s) / s);

  if (cons_P > 1e-3 || cons_s > 1e-3) {
    cout << "TD consistency not verified in hs call" << endl;
  }
}

su2double CNonPolytropicPengRobinson::ComputeIntegralCp0_DT(su2double T){
  su2double int_cp0 = 0;

  for(size_t i_cp=0; i_cp < coeffs_.size(); i_cp++){
    int_cp0 += coeffs_.at(i_cp) * pow(T, (i_cp + 1.0)) / (i_cp + 1.0);
  }
  return int_cp0;
}

su2double CNonPolytropicPengRobinson::ComputeIntegralCp0_T_DT(su2double T){
  su2double int_cp0_T = 0;

  for(size_t i_cp=0; i_cp < coeffs_.size(); i_cp++){

    if (i_cp == 0) {
      int_cp0_T += coeffs_.at(i_cp) * log(T);
    } else {
      int_cp0_T += coeffs_.at(i_cp) * pow(T, i_cp) / (i_cp);
    }
  }
  return int_cp0_T;
}

su2double CNonPolytropicPengRobinson::GetCv0(){
  su2double cv0, cp0=0;

  for(size_t i_cp=0; i_cp < coeffs_.size(); i_cp++){
    cp0 += coeffs_.at(i_cp) * pow(Temperature, i_cp);
  }

  cv0 = cp0 - Gas_Constant;

  return cv0;
}