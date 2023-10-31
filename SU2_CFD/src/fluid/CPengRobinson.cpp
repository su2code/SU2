/*!
 * \file CPengRobinson.cpp
 * \brief Source of the Peng-Robinson model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CPengRobinson.hpp"

CPengRobinson::CPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w)
    : CIdealGas(gamma, R) {
  a = 0.45724 * Gas_Constant * Gas_Constant * Tstar * Tstar / Pstar;
  b = 0.0778 * Gas_Constant * Tstar / Pstar;
  TstarCrit = Tstar;
  Zed = 1.0;

  if (w <= 0.49)
    k = 0.37464 + 1.54226 * w - 0.26992 * w * w;
  else
    k = 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w * w * w;
}

su2double CPengRobinson::alpha2(su2double T) const {
  return (1 + k * (1 - sqrt(T / TstarCrit))) * (1 + k * (1 - sqrt(T / TstarCrit)));
}

su2double CPengRobinson::T_v_h(su2double v, su2double h) {
  su2double fv, A, B, C, T, d, atanh;
  su2double sqrt2 = sqrt(2.0);

  d = (v * v + 2 * b * v - b * b);

  atanh = (log(1.0 + (b * sqrt2 / (v + b))) - log(1.0 - (b * sqrt2 / (v + b)))) / 2.0;

  fv = atanh;

  A = Gas_Constant * (1 / Gamma_Minus_One + v / (v - b)) - a * v * k * k / (TstarCrit * d);
  B = a * k * (k + 1) / sqrt(TstarCrit) * (fv / (b * sqrt2) + 2 * v / d);
  C = h + a * (1 + k) * (1 + k) * (fv / (b * sqrt2) + v / d);

  T = (-B + sqrt(B * B + 4 * A * C)) / (2 * A);  /// Only positive root considered

  return T * T;
}

su2double CPengRobinson::T_P_rho(su2double P, su2double rho) {
  su2double A, B, C, T, vb1, vb2;
  vb1 = (1 / rho - b);
  vb2 = (1 / rho / rho + 2 * b / rho - b * b);

  A = Gas_Constant / vb1 - a * k * k / TstarCrit / vb2;

  B = 2 * a * k * (k + 1) / sqrt(TstarCrit) / vb2;

  C = -P - a * (1 + k) * (1 + k) / vb2;

  T = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
  T *= T;
  return T;
}

void CPengRobinson::SetTDState_rhoe(su2double rho, su2double e) {
  su2double DpDd_T, DpDT_d, DeDd_T, Cv;
  su2double A, B, C, sqrt2, fv, a2T, rho2, atanh;

  Density = rho;
  StaticEnergy = e;

  AD::StartPreacc();
  AD::SetPreaccIn(rho);
  AD::SetPreaccIn(e);

  rho2 = rho * rho;
  sqrt2 = sqrt(2.0);

  atanh = (log(1.0 + (rho * b * sqrt2 / (1 + rho * b))) - log(1.0 - (rho * b * sqrt2 / (1 + rho * b)))) / 2.0;

  fv = atanh;

  A = Gas_Constant / Gamma_Minus_One;
  B = a * k * (k + 1) * fv / (b * sqrt2 * sqrt(TstarCrit));
  C = a * (k + 1) * (k + 1) * fv / (b * sqrt2) + e;

  Temperature = (-B + sqrt(B * B + 4 * A * C)) / (2 * A);  /// Only positive root considered
  Temperature *= Temperature;

  a2T = alpha2(Temperature);

  A = (1 / rho2 + 2 * b / rho - b * b);
  B = 1 / rho - b;

  Pressure = Temperature * Gas_Constant / B - a * a2T / A;

  Entropy = Gas_Constant / Gamma_Minus_One * log(Temperature) + Gas_Constant * log(B) -
            a * sqrt(a2T) * k * fv / (b * sqrt2 * sqrt(Temperature * TstarCrit));

  DpDd_T = (Temperature * Gas_Constant / (B * B) - 2 * a * a2T * (1 / rho + b) / (A * A)) / (rho2);

  DpDT_d = Gas_Constant / B + a * k / A * sqrt(a2T / (Temperature * TstarCrit));

  Cv = Gas_Constant / Gamma_Minus_One + (a * k * (k + 1) * fv) / (2 * b * sqrt(2 * Temperature * TstarCrit));

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

void CPengRobinson::SetTDState_PT(su2double P, su2double T) {
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

  e = T * Gas_Constant / Gamma_Minus_One - a * (k + 1) * sqrt(alpha2(T)) * fv / (b * sqrt2);

  AD::SetPreaccOut(rho);
  AD::SetPreaccOut(e);
  AD::EndPreacc();

  SetTDState_rhoe(rho, e);
}

void CPengRobinson::SetTDState_Prho(su2double P, su2double rho) {
  SetEnergy_Prho(P, rho);

  SetTDState_rhoe(rho, StaticEnergy);
}

void CPengRobinson::SetTDState_hs(su2double h, su2double s) {
  su2double T, fv, sqrt2 = sqrt(2.0), A;
  su2double f, v, atanh;
  su2double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
  su2double toll = 1e-9, FACTOR = 0.2;
  su2double cons_s, cons_h;
  unsigned short countrtb = 0, NTRY = 100, ITMAX = 100;

  A = Gas_Constant / Gamma_Minus_One;
  T = h * Gamma_Minus_One / Gas_Constant / Gamma;
  v = exp(-1 / Gamma_Minus_One * log(T) + s / Gas_Constant);

  x1 = 0.2 * v;
  x2 = 0.35 * v;

  T = T_v_h(x1, h);

  atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
  fv = atanh;

  fx1 = A * log(T) + Gas_Constant * log(x1 - b) - a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
  T = T_v_h(x2, h);

  atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
  fv = atanh;

  fx2 = A * log(T) + Gas_Constant * log(x2 - b) - a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;

  // zbrac algorithm NR

  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        T = T_v_h(x1, h);
        atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
        fv = atanh;
        fx1 = A * log(T) + Gas_Constant * log(x1 - b) -
              a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
      } else {
        x2 += FACTOR * (x2 - x1);
        T = T_v_h(x2, h);
        atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
        fv = atanh;
        fx2 = A * log(T) + Gas_Constant * log(x2 - b) -
              a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
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
    fmid = A * log(T) + Gas_Constant * log(xmid - b) -
           a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;

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

void CPengRobinson::SetEnergy_Prho(su2double P, su2double rho) {
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

  StaticEnergy = T * Gas_Constant / Gamma_Minus_One - ad;

  AD::SetPreaccOut(StaticEnergy);
  AD::EndPreacc();
}

void CPengRobinson::SetTDState_rhoT(su2double rho, su2double T) {
  su2double fv, e, atanh;

  atanh = (log(1.0 + (rho * b * sqrt(2.0) / (1 + rho * b))) - log(1.0 - (rho * b * sqrt(2.0) / (1 + rho * b)))) / 2.0;
  fv = atanh;
  e = T * Gas_Constant / Gamma_Minus_One - a * (k + 1) * sqrt(alpha2(T)) / (b * sqrt(2.0)) * fv;
  SetTDState_rhoe(rho, e);
}

void CPengRobinson::SetTDState_Ps(su2double P, su2double s) {
  su2double T, rho, v, cons_P, cons_s, fv, A, atanh;
  su2double x1, x2, fx1, fx2, f, fmid, rtb, dx, xmid, sqrt2 = sqrt(2.0);
  su2double toll = 1e-5, FACTOR = 0.2;
  unsigned short count = 0, NTRY = 100, ITMAX = 100;

  A = Gas_Constant / Gamma_Minus_One;
  T = exp(Gamma_Minus_One / Gamma * (s / Gas_Constant + log(P) - log(Gas_Constant)));
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

  fx1 = A * log(T) + Gas_Constant * log(x1 - b) - a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
  T = T_P_rho(P, 1.0 / x2);

  atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
  fv = atanh;

  fx2 = A * log(T) + Gas_Constant * log(x2 - b) - a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;

  // zbrac algorithm NR

  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        T = T_P_rho(P, 1.0 / x1);

        atanh = (log(1.0 + (b * sqrt2 / (x1 + b))) - log(1.0 - (b * sqrt2 / (x1 + b)))) / 2.0;
        fv = atanh;

        fx1 = A * log(T) + Gas_Constant * log(x1 - b) -
              a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
      } else {
        T = T_P_rho(P, 1.0 / x2);

        atanh = (log(1.0 + (b * sqrt2 / (x2 + b))) - log(1.0 - (b * sqrt2 / (x2 + b)))) / 2.0;
        fv = atanh;

        fx2 = A * log(T) + Gas_Constant * log(x2 - b) -
              a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
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

    fmid = A * log(T) + Gas_Constant * log(xmid - b) -
           a * sqrt(alpha2(T)) * k * fv / (b * sqrt2 * sqrt(T * TstarCrit)) - s;
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

void CPengRobinson::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  su2double dPdT_rho, dPdrho_T, dPds_rho, der1_alpha;

  SetTDState_Prho(P, rho);

  der1_alpha = -k / (2 * sqrt(TstarCrit * Temperature));
  dPdT_rho = Gas_Constant * rho / (1.0 - rho * b) -
             2 * a * rho * rho * sqrt(alpha2(Temperature)) * der1_alpha / (1 + 2 * b * rho - b * b * rho * rho);
  dPdrho_T = Gas_Constant * Temperature / (1.0 - rho * b) / (1.0 - rho * b) -
             2.0 * rho * a * alpha2(Temperature) * (1.0 + b * rho) / (1 + 2 * b * rho - b * b * rho * rho) /
                 (1 + 2 * b * rho - b * b * rho * rho);

  dhdrho_P = -dPdrho_e / dPde_rho - P / rho / rho;
  dhdP_rho = 1.0 / dPde_rho + 1.0 / rho;
  dPds_rho = rho * rho * (SoundSpeed2 - dPdrho_T) / dPdT_rho;
  dsdP_rho = 1.0 / dPds_rho;
  dsdrho_P = -SoundSpeed2 / dPds_rho;
}
