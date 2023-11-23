/*!
 * \file C1DInterpolation.cpp
 * \brief Classes for 1D interpolation.
 * \author Aman Baig, P. Gomes
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

#include "../../include/toolboxes/C1DInterpolation.hpp"
#include "../../include/linear_algebra/blas_structure.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

su2double CAkimaInterpolation::EvaluateSpline(su2double Point_Interp) const {
  const auto i = lower_bound(Point_Interp);

  if (i >= x.size() - 1) return (Point_Interp <= x[0]) ? y.front() : y.back();

  const su2double h = Point_Interp - x[i];

  return y[i] + h * (b[i] + h * (c[i] + h * d[i]));
}

su2double CLinearInterpolation::EvaluateSpline(su2double Point_Interp) const {
  const auto i = lower_bound(Point_Interp);

  if (i >= x.size() - 1) return (Point_Interp <= x[0]) ? y.front() : y.back();

  return y[i] + (Point_Interp - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

void CCubicSpline::SetSpline(const vector<su2double>& X, const vector<su2double>& Data) {
  C1DInterpolation::SetSpline(X, Data);

  const int N = x.size();

  /*--- Alias the vectors of coefficients to build the tridiagonal system. ---*/

  auto& lower = b;
  b.resize(N);
  auto& main = c;
  c.resize(N);
  auto& upper = d;
  d.resize(N);
  vector<su2double> rhs(N);

  /*--- Main part of the tridiagonal system. ---*/

  for (int i = 1; i < N - 1; i++) {
    lower[i] = x[i] - x[i - 1];
    main[i] = (x[i + 1] - x[i - 1]) * 2;
    upper[i] = x[i + 1] - x[i];
    rhs[i] = 6 * ((y[i + 1] - y[i]) / upper[i] - (y[i] - y[i - 1]) / lower[i]);
  }

  /*--- Start condition. ---*/

  if (startDer == SECOND) {
    main[0] = 1.0;
    upper[0] = 0;
    rhs[0] = startVal;
  } else {  // FIRST
    main[0] = 2 * lower[1];
    upper[0] = lower[1];
    rhs[0] = 6 * ((y[1] - y[0]) / lower[1] - startVal);
  }

  /*--- End condition. ---*/

  if (endDer == SECOND) {
    main[N - 1] = 1.0;
    lower[N - 1] = 0;
    rhs[N - 1] = endVal;
  } else {  // FIRST
    main[N - 1] = 2 * upper[N - 2];
    lower[N - 1] = upper[N - 2];
    rhs[N - 1] = 6 * (endVal - (y[N - 1] - y[N - 2]) / upper[N - 2]);
  }

  /*--- Solve system for 2nd derivative at the knots. ---*/

  CBlasStructure::tdma(lower, main, upper, rhs);
  const auto& d2y = rhs;

  /*--- Compute the polynomial coefficients. ---*/

  for (int i = 0; i < N - 1; i++) {
    su2double dx = x[i + 1] - x[i];
    c[i] = 0.5 * d2y[i];
    d[i] = (d2y[i + 1] - d2y[i]) / (6 * dx);
    b[i] = (y[i + 1] - y[i]) / dx - (c[i] + d[i] * dx) * dx;
  }
}

void CAkimaInterpolation::SetSpline(const vector<su2double>& X, const vector<su2double>& Data) {
  C1DInterpolation::SetSpline(X, Data);

  const int n = X.size();
  vector<su2double> h(n - 1);
  vector<su2double> p(n - 1);

  /*---calculating finite differences (h) and gradients (p) ---*/
  for (int i = 0; i < n - 1; i++) {
    h[i] = X[i + 1] - X[i];
    p[i] = (Data[i + 1] - Data[i]) / h[i];
  }

  /*--- b,c,d are the akima spline's cofficient for the cubic equation ---*/
  b.resize(n);
  c.resize(n);
  d.resize(n);

  b[0] = p[0];
  b[1] = (p[0] + p[1]) / 2;
  b[n - 1] = p[n - 2];
  b[n - 2] = (p[n - 2] + p[n - 3]) / 2;

  for (int i = 2; i < n - 2; i++) {
    su2double w1 = fabs(p[i + 1] - p[i]), w2 = fabs(p[i - 1] - p[i - 2]);
    if (w1 + w2 < 0.0001) {
      b[i] = (p[i - 1] + p[i]) / 2;
    } else {
      b[i] = (w1 * p[i - 1] + w2 * p[i]) / (w1 + w2);
    }
  }

  for (int i = 0; i < n - 1; i++) {
    c[i] = (3 * p[i] - 2 * b[i] - b[i + 1]) / h[i];
    d[i] = (b[i + 1] + b[i] - 2 * p[i]) / h[i] / h[i];
  }
  b.back() = c.back() = 0;
}

vector<su2double> CorrectedInletValues(const vector<su2double>& Inlet_Interpolated, su2double Theta,
                                       unsigned short nDim, const su2double* Coord, unsigned short nVar_Turb,
                                       INLET_INTERP_TYPE Interpolation_Type) {
  unsigned short size_columns = Inlet_Interpolated.size() + nDim;
  vector<su2double> Inlet_Values(size_columns);
  su2double unit_r, unit_Theta, unit_m, Alpha, Phi;

  /*---For x,y,z,T,P columns---*/
  for (int i = 0; i < nDim; i++) Inlet_Values[i] = Coord[i];

  for (int i = nDim; i < nDim + 2; i++) Inlet_Values[i] = Inlet_Interpolated[i - 2];

  /*---For turbulence variables columns---*/
  if (nVar_Turb == 1)
    Inlet_Values[nDim + 5] = Inlet_Interpolated[5];
  else if (nVar_Turb == 2) {
    Inlet_Values[nDim + 5] = Inlet_Interpolated[5];
    Inlet_Values[nDim + 6] = Inlet_Interpolated[6];
  }

  /*--- Correct for Interpolation Type now ---*/
  switch (Interpolation_Type) {
    case (INLET_INTERP_TYPE::VR_VTHETA):
      unit_r = Inlet_Interpolated[nDim];
      unit_Theta = Inlet_Interpolated[nDim + 1];
      break;

    case (INLET_INTERP_TYPE::ALPHA_PHI):
      Alpha = Inlet_Interpolated[nDim] * PI_NUMBER / 180;
      Phi = Inlet_Interpolated[nDim + 1] * PI_NUMBER / 180;
      unit_m = sqrt(1 / (1 + pow(tan(Alpha), 2)));
      unit_Theta = tan(Alpha) * unit_m;
      unit_r = unit_m * sin(Phi);
      break;

    default:
      unit_Theta = unit_r = 0.0;
      break;
  }

  /*--- Converting from cylindrical to cartesian unit vectors ---*/
  Inlet_Values[nDim + 2] = unit_r * cos(Theta) - unit_Theta * sin(Theta);
  Inlet_Values[nDim + 3] = unit_r * sin(Theta) + unit_Theta * cos(Theta);
  Inlet_Values[nDim + 4] = sqrt(1 - pow(unit_r, 2) - pow(unit_Theta, 2));

  return Inlet_Values;
}

void PrintInletInterpolatedData(const vector<su2double>& Inlet_Data_Interpolated, const string& Marker,
                                unsigned long nVertex, unsigned short nDim, unsigned short nColumns) {
  ofstream myfile;
  myfile.precision(16);
  myfile.open("Interpolated_Data_" + Marker + ".dat", ios_base::out);

  if (myfile.is_open()) {
    for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
      for (unsigned short iVar = 0; iVar < nColumns; iVar++) {
        myfile << Inlet_Data_Interpolated[iVertex * nColumns + iVar] << "\t";
      }
      myfile << endl;
    }
    myfile.close();
  } else
    cout << "file cannot be opened" << endl;
}
