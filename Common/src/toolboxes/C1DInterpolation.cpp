/*!
 * \file C1DInterpolation.cpp
 * \brief Inlet_interpolation_functions
 * \author Aman Baig
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

su2double CAkimaInterpolation::EvaluateSpline(su2double Point_Interp){

  Point_Match = false;

  for (int i=0; i<n-1; i++)
    if(Point_Interp>=x[i] && Point_Interp<=x[i+1])
      Point_Match = true;

  if (Point_Match == false){
    cout<<"WARNING: Extrapolating data for radius: "<<Point_Interp<<endl;
    Point_Match = true;
  }

  int i=0, j=n-1;

  while (j-i>1){
    int m = (i+j) / 2;
    if (Point_Interp>x[m]) {i = m;}
    else {j = m;}
  }

  su2double h=Point_Interp-x[i];

  return y[i]+h*(b[i]+h*(c[i]+h*d[i]));
}

su2double CLinearInterpolation::EvaluateSpline(su2double Point_Interp){

  Point_Match = false;

  for (int i=0; i<n-1; i++){
    if(Point_Interp>=x[i] && Point_Interp<=x[i+1]){
      Point_Match = true;
      return (Point_Interp-x[i])*dydx[i]+y[i];
    }
  }
  return 0;
}

void CLinearInterpolation::SetSpline(const vector<su2double> &X, const vector<su2double> &Data){

  n = X.size();
  su2double h;
  x.resize(n);
  y.resize(n);
  dydx.resize(n-1);

  x = X;
  y = Data;

  for (int i=0; i<n-1;i ++){
    h = x[i+1]-x[i];
    dydx[i] = (y[i+1]-y[i])/h;
  }
}

void CAkimaInterpolation::SetSpline(const vector<su2double> &X,const vector<su2double> &Data){

  n = X.size();
  vector<su2double> h (n-1);
  vector<su2double> p (n-1);

  /*---calculating finite differences (h) and gradients (p) ---*/
  for (int i=0; i<n-1;i++){
      h[i]=X[i+1]-X[i];
      p[i]=(Data[i+1]-Data[i])/h[i];
      }

  /*---b,c,d are the akima spline's cofficient for the cubic equation---*/
  x.resize(n);
  y.resize(n);
  b.resize(n);
  c.resize(n-1);
  d.resize(n-1);

  x=X;
  y=Data;

  b[0] = p[0] ;
  b[1] =(p[0]+ p[1])/2 ;
  b[n-1]=p[n-2];
  b[n-2]=(p[n-2]+p[n-3])/2;

  for (int i = 2; i < n-2; i ++){
    su2double w1=fabs(p[i+1]-p[i]) , w2=fabs(p[i-1]-p[i-2]);
    if (w1+w2<0.0001) {
      b[i] = (p[i-1]+p[i])/2 ;
    }
    else{
      b[i] = (w1*p[i-1]+w2*p[i])/(w1+w2);
    }
  }

  for (int i = 0; i < n-1; i ++){
    c[i] = (3*p[i]-2*b[i]-b[i+1])/h[i] ;
    d[i] = (b[i+1]+b[i]-2*p[i])/h[i]/h[i];
  }
}

vector<su2double> CorrectedInletValues(const vector<su2double> &Inlet_Interpolated ,
                                       su2double Theta ,
                                       unsigned short nDim,
                                       const su2double *Coord,
                                       unsigned short nVar_Turb,
                                       ENUM_INLET_INTERPOLATIONTYPE Interpolation_Type){

  unsigned short size_columns=Inlet_Interpolated.size()+nDim;
  vector<su2double> Inlet_Values(size_columns);
  su2double unit_r, unit_Theta, unit_m, Alpha, Phi;

  /*---For x,y,z,T,P columns---*/
  for (int i=0;i<nDim;i++)
    Inlet_Values[i] = Coord[i];

  for (int i=nDim;i<nDim+2;i++)
    Inlet_Values[i] = Inlet_Interpolated[i-2];

  /*---For turbulence variables columns---*/
  if (nVar_Turb == 1)
    Inlet_Values[nDim+5] = Inlet_Interpolated[5];
  else if (nVar_Turb == 2){
    Inlet_Values[nDim+5] = Inlet_Interpolated[5];
    Inlet_Values[nDim+6] = Inlet_Interpolated[6];
  }

  /*--- Correct for Interpolation Type now ---*/
  switch(Interpolation_Type){
    case(VR_VTHETA):
      unit_r = Inlet_Interpolated[nDim];
      unit_Theta = Inlet_Interpolated[nDim+1];
      break;

    case(ALPHA_PHI):
      Alpha = Inlet_Interpolated[nDim]*PI_NUMBER/180;
      Phi = Inlet_Interpolated[nDim+1]*PI_NUMBER/180;
      unit_m = sqrt(1/(1+pow(tan(Alpha),2)));
      unit_Theta = tan(Alpha)*unit_m;
      unit_r=unit_m*sin(Phi);
      break;

    default:
      unit_Theta = unit_r = 0.0;
      break;
  }

  /*--- Converting from cylindrical to cartesian unit vectors ---*/
  Inlet_Values[nDim+2] = unit_r*cos(Theta) - unit_Theta*sin(Theta);
  Inlet_Values[nDim+3] = unit_r*sin(Theta) + unit_Theta*cos(Theta);
  Inlet_Values[nDim+4] = sqrt(1-pow(unit_r,2)- pow(unit_Theta,2));

  return Inlet_Values;
}

void PrintInletInterpolatedData(const vector<su2double>& Inlet_Data_Interpolated, string Marker,
                                unsigned long nVertex, unsigned short nDim, unsigned short nColumns){

  ofstream myfile;
  myfile.precision(16);
  myfile.open("Interpolated_Data_"+Marker+".dat",ios_base::out);

  if (myfile.is_open()){
    for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
      for (unsigned short iVar=0; iVar < nColumns; iVar++){
        myfile<<Inlet_Data_Interpolated[iVertex*nColumns+iVar]<<"\t";
      }
      myfile<<endl;
    }
    myfile.close();
  }
  else
    cout<<"file cannot be opened"<<endl;
}
