/*!
 * \file CVanDerWaalsGas.cpp
 * \brief Source of the Polytropic Van der Waals model.
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

#include "../../include/fluid/CVanDerWaalsGas.hpp"

CVanDerWaalsGas::CVanDerWaalsGas(su2double gamma, su2double R, su2double Pstar, su2double Tstar) : CIdealGas(gamma, R) {
  a = 27.0 / 64.0 * Gas_Constant * Gas_Constant * Tstar * Tstar / Pstar;
  b = 1.0 / 8.0 * Gas_Constant * Tstar / Pstar;
  Zed = 1.0;
}

void CVanDerWaalsGas::SetTDState_rhoe(su2double rho, su2double e) {
  Density = rho;
  StaticEnergy = e;

  Pressure = Gamma_Minus_One * Density / (1.0 - Density * b) * (StaticEnergy + Density * a) - a * Density * Density;
  Temperature = (Pressure + Density * Density * a) * ((1 - Density * b) / (Density * Gas_Constant));
  Entropy = Gas_Constant * (log(Temperature) / Gamma_Minus_One + log(1 / Density - b));

  dPde_rho = Density * Gamma_Minus_One / (1.0 - Density * b);
  dPdrho_e = Gamma_Minus_One / (1.0 - Density * b) *
                 ((StaticEnergy + 2 * Density * a) + Density * b * (StaticEnergy + Density * a) / (1.0 - Density * b)) -
             2 * Density * a;
  dTdrho_e = Gamma_Minus_One / Gas_Constant * a;
  dTde_rho = Gamma_Minus_One / Gas_Constant;

  SoundSpeed2 = dPdrho_e + Pressure / (Density * Density) * dPde_rho;

  Zed = Pressure / (Gas_Constant * Temperature * Density);
}

void CVanDerWaalsGas::SetTDState_PT(su2double P, su2double T) {
  su2double toll = 1e-5;
  unsigned short nmax = 20, count = 0;
  su2double A, B, Z, DZ = 1.0, F, F1;
  A = a * P / (T * Gas_Constant) / (T * Gas_Constant);
  B = b * P / (T * Gas_Constant);

  if (Zed > 0.1)
    Z = min(Zed, 0.99);
  else
    Z = 0.99;

  do {
    F = Z * Z * Z - Z * Z * (B + 1.0) + Z * A - A * B;
    F1 = 3 * Z * Z - 2 * Z * (B + 1.0) + A;
    DZ = F / F1;
    Z -= 0.7 * DZ;
    count++;
  } while (abs(DZ) > toll && count < nmax);

  if (count == nmax) {
    cout << "Warning Newton-Raphson exceed number of max iteration in PT" << endl;
    cout << "Compressibility factor  " << Z << " would be substituted with " << Zed << endl;
  }

  // check if the solution is physical otherwise uses previous point solution
  if (Z <= 1.01 && Z >= 0.05 && count < nmax) Zed = Z;

  Density = P / (Zed * Gas_Constant * T);

  su2double e = T * Gas_Constant / Gamma_Minus_One - a * Density;
  SetTDState_rhoe(Density, e);
}

void CVanDerWaalsGas::SetTDState_Prho(su2double P, su2double rho) {
  SetEnergy_Prho(P, rho);
  SetTDState_rhoe(rho, StaticEnergy);
}

void CVanDerWaalsGas::SetTDState_hs(su2double h, su2double s) {
  su2double v, T, rho, f, fmid, rtb;
  su2double x1, x2, xmid, dx, fx1, fx2;
  su2double toll = 1e-5, FACTOR = 0.2;
  unsigned short count = 0, NTRY = 100, ITMAX = 100;
  su2double cons_s, cons_h;

  T = 1.0 * h * Gamma_Minus_One / Gas_Constant / Gamma;
  v = exp(-1 / Gamma_Minus_One * log(T) + s / Gas_Constant);
  x1 = 0.2 * v;
  x2 = 0.35 * v;
  fx1 = log(x1 - b) - s / Gas_Constant +
        log((h + 2 * a / x1) / Gas_Constant / (1 / Gamma_Minus_One + x1 / (x1 - b))) / Gamma_Minus_One;
  fx2 = log(x2 - b) - s / Gas_Constant +
        log((h + 2 * a / x2) / Gas_Constant / (1 / Gamma_Minus_One + x2 / (x2 - b))) / Gamma_Minus_One;

  // zbrac algorithm NR
  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        fx1 = log(x1 - b) - s / Gas_Constant +
              log((h + 2 * a / x1) / Gas_Constant / (1 / Gamma_Minus_One + x1 / (x1 - b))) / Gamma_Minus_One;
      } else {
        x2 += FACTOR * (x2 - x1);
        fx2 = log(x2 - b) - s / Gas_Constant +
              log((h + 2 * a / x2) / Gas_Constant / (1 / Gamma_Minus_One + x2 / (x2 - b))) / Gamma_Minus_One;
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
    fmid = log(xmid - b) - s / Gas_Constant +
           log((h + 2 * a / xmid) / Gas_Constant / (1 / Gamma_Minus_One + xmid / (xmid - b))) / Gamma_Minus_One;
    if (fmid <= 0.0) rtb = xmid;
    count++;
  } while (abs(fmid) > toll && count < ITMAX);

  v = xmid;
  if (count == ITMAX) {
    cout << "Too many bisections in rtbis" << endl;
  }

  rho = 1 / v;
  T = (h + 2 * a / v) / Gas_Constant / (1 / Gamma_Minus_One + v / (v - b));
  SetTDState_rhoT(rho, T);

  cons_h = abs(((StaticEnergy + Pressure / Density) - h) / h);
  cons_s = abs((Entropy - s) / s);

  if (cons_h > 1e-3 || cons_s > 1e-3) {
    cout << "TD consistency not verified in hs call" << endl;
  }
}

void CVanDerWaalsGas::SetEnergy_Prho(su2double P, su2double rho) {
  su2double T = (P + rho * rho * a) * (1 - rho * b) / (rho * Gas_Constant);
  StaticEnergy = T * Gas_Constant / Gamma_Minus_One - rho * a;
}

void CVanDerWaalsGas::SetTDState_rhoT(su2double rho, su2double T) {
  su2double e = T * Gas_Constant / Gamma_Minus_One - a * rho;
  SetTDState_rhoe(rho, e);
}

void CVanDerWaalsGas::SetTDState_Ps(su2double P, su2double s) {
  su2double T, rho, cons_P, cons_s;
  su2double x1, x2, fx1, fx2, f, fmid, T1, T2, rtb, dx, xmid;
  su2double toll = 1e-5, FACTOR = 0.2;
  unsigned short count = 0, NTRY = 100, ITMAX = 100;

  T = exp(Gamma_Minus_One / Gamma * (s / Gas_Constant + log(P) - log(Gas_Constant)));
  rho = P / (T * Gas_Constant);

  if (Zed < 0.9999) {
    x1 = rho;
    x2 = rho / Zed;

  } else {
    x1 = rho;
    x2 = rho / 0.5;
  }
  T1 = (P + x1 * x1 * a) * ((1 - x1 * b) / (x1 * Gas_Constant));
  fx1 = Gas_Constant * (log(T1) / Gamma_Minus_One + log(1 / x1 - b)) - s;
  T2 = (P + x2 * x2 * a) * ((1 - x2 * b) / (x2 * Gas_Constant));
  fx2 = Gas_Constant * (log(T2) / Gamma_Minus_One + log(1 / x2 - b)) - s;

  // zbrac algorithm NR
  for (int j = 1; j <= NTRY; j++) {
    if (fx1 * fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR * (x1 - x2);
        T1 = (P + x1 * x1 * a) * ((1 - x1 * b) / (x1 * Gas_Constant));
        fx1 = Gas_Constant * (log(T1) / Gamma_Minus_One + log(1 / x1 - b)) - s;
      } else {
        x2 += FACTOR * (x2 - x1);
        T2 = (P + x2 * x2 * a) * ((1 - x2 * b) / (x2 * Gas_Constant));
        fx2 = Gas_Constant * (log(T2) / Gamma_Minus_One + log(1 / x2 - b)) - s;
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
    T = (P + xmid * xmid * a) * ((1 - xmid * b) / (xmid * Gas_Constant));
    fmid = Gas_Constant * (log(T) / Gamma_Minus_One + log(1 / xmid - b)) - s;
    if (fmid <= 0.0) rtb = xmid;
    count++;
  } while (abs(fmid) > toll && count < ITMAX);

  if (count == ITMAX) {
    cout << "Too many bisections in rtbis" << endl;
  }

  rho = xmid;
  T = (P + rho * rho * a) * ((1 - rho * b) / (rho * Gas_Constant));
  SetTDState_rhoT(rho, T);

  cons_P = abs((Pressure - P) / P);
  cons_s = abs((Entropy - s) / s);

  if (cons_P > 1e-3 || cons_s > 1e-3) {
    cout << "TD consistency not verified in hs call" << endl;
  }
}

void CVanDerWaalsGas::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  su2double dPdT_rho, dPdrho_T, dPds_rho;

  SetTDState_Prho(P, rho);

  dPdT_rho = Gas_Constant * rho / (1.0 - rho * b);
  dPdrho_T = Gas_Constant * Temperature / (1.0 - rho * b) / (1.0 - rho * b) - 2.0 * rho * a;

  dhdrho_P = -dPdrho_e / dPde_rho - P / rho / rho;
  dhdP_rho = 1.0 / dPde_rho + 1.0 / rho;
  dPds_rho = rho * rho * (SoundSpeed2 - dPdrho_T) / dPdT_rho;
  dsdP_rho = 1.0 / dPds_rho;
  dsdrho_P = -SoundSpeed2 / dPds_rho;
}
