/*!
 * \file inlet_interpolation_functions.cpp
 * \brief Inlet_interpolation_functions
 * \author Aman Baig
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
 

#include "../../include/toolboxes/CInletInterpolation.hpp"

void CAkimaSpline::akima_spline_set (CAkimaSpline *s, vector<su2double> &x,vector<su2double> &y){
int n = x.size();
vector<su2double> h (n-1);
vector<su2double> p (n-1);

for (int i=0; i<n-1;i ++){h[i]=x[i+1]-x[i];}
for (int i=0; i<n-1;i++) p[i]=(y[i+1]-y[i]) / h [i] ;

s->x.resize(n);
s->y.resize(n);
s->b.resize(n);
s->c.resize(n-1);
s->d.resize(n-1);

s->n=n; for (int i=0; i<n ; i++){s->x[i]=x[i]; s->y[i]=y[i];}
s->b[0] = p[0] ; s->b[1] =(p[0]+ p[1])/2 ;

s->b[n-1]=p[n-2]; s->b [n-2]=(p[n-2]+p[n-3])/2;

for (int i =2; i<n-2; i ++){
su2double w1=fabs(p[i+1]-p[i]) , w2=fabs(p[i-1]-p[i-2]);
if (w1+w2==0) s->b[i] = (p[i-1]+p[i])/2 ;
else s->b [i] = (w1*p[i-1]+w2*p[i])/(w1+w2) ;
}
for (int i =0; i<n-1; i ++){
s->c [i]=(3*p [i]-2*s->b [ i ]-s->b[i+1])/h[i] ;
s->d [i]=(s->b[i+1]+s->b [ i ]-2*p[i])/h[i]/h[i];
}
}

su2double CAkimaSpline::akima_spline_eval(CAkimaSpline *s, su2double Point_Interp){
int i =0, j=s->n-1;
while ( j-i >1){ int m=( i+j ) / 2 ; if ( Point_Interp>s->x[m] ) i=m; else j=m; }
su2double h=Point_Interp-s->x[i] ;
return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i])) ;
}

void CLinearSpline::Linear_set(CLinearSpline *s,vector<su2double> &x, vector<su2double> &y){
int n = x.size();
vector<su2double> h (n-1);
s->x.resize(n);
s->y.resize(n);
s->dydx.resize(n-1);
s->x = x;
s->y = y;
for (int i=0; i<n-1;i ++){h[i]=x[i+1]-x[i];}
for (int i=0; i<n-1;i++) {s->dydx[i]=(y[i+1]-y[i]) / h [i];}
}

su2double CLinearSpline::Linear_eval(CLinearSpline *s, su2double Point_Interp){
    int size = s->x.size();

    for (int i=0;i<size-1;i++){
        if(Point_Interp>=s->x[i] && Point_Interp<=s->x[i+1]){
            s->Point_Match = true;
            return (Point_Interp-s->x[i])*dydx[i]+s->y[i];}
    }
}

vector<su2double> CInletInterpolation::CorrectedInletValues(vector<su2double> Inlet_Interpolated , 
                                                    su2double Theta ,
                                                    unsigned short nDim, 
                                                    su2double *Coord, 
                                                    unsigned short nVar_Turb,
                                                    string Interpolation_Type){
su2double size_columns=Inlet_Interpolated.size();
vector<su2double> Inlet_Values (size_columns+nDim);
su2double unit_r, unit_Theta, unit_m, Alpha, Phi;

/*---For x,y,z,T,P columns---*/
for (int i=0;i<nDim+2;i++){
    if (i<nDim)
        Inlet_Values[i] = Coord[i];
    else
        Inlet_Values[i] = Inlet_Interpolated[nDim-2];
}

/*---For turbulence variables columns---*/
if (nVar_Turb == 1)
    Inlet_Values[nDim+5] = Inlet_Values[5];
if (nVar_Turb == 2)
    Inlet_Values[nDim+6] = Inlet_Values[6];

/*--- Correct for Interpolation Type now ---*/
if(Interpolation_Type == "VRVTHETA"){
    unit_r = Inlet_Interpolated[nDim+2];
    unit_Theta = Inlet_Interpolated[nDim+3];
}
else if(Interpolation_Type == "ALPHA_PHI"){
    Alpha = Inlet_Interpolated[nDim+2];
    Phi = Inlet_Interpolated[nDim+3];
    unit_m = sqrt(1/(1+pow(tan(Alpha),2)));

    unit_Theta = tan(Alpha)*unit_m;
    unit_r=unit_m*sin(Phi);
}

Inlet_Values[nDim+2] = unit_r*cos(Theta) - unit_Theta*sin(Theta); //for ix
Inlet_Values[nDim+3] = unit_r*sin(Theta) + unit_Theta*cos(Theta); //for iy
Inlet_Values[nDim+4] = sqrt(1-pow(unit_r,2)- pow(unit_Theta,2));  //for iz

return Inlet_Values;

}

void CInletInterpolation::PrintInletInterpolatedData(vector<su2double> Inlet_Values, string Marker, unsigned long nVertex, unsigned short nDim){
ofstream myfile;
myfile.precision(16);
myfile.open("Interpolated_Data_"+Marker+".dat",ios_base::out);

if (myfile.is_open()){  
  for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
            for  (unsigned long iVar=0; iVar < Inlet_Values.size(); iVar++){
                myfile<<Inlet_Values[iVertex*Inlet_Values.size()+iVar]<<"\t";
              }
            myfile<<endl;
          }
  myfile.close();
}
  else
    cout<<"file cannot be opened"<<endl;
}
