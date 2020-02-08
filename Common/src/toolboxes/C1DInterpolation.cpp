/*!
 * \file inlet_interpolation_functions.cpp
 * \brief Inlet_interpolation_functions
 * \author Aman Baig
 * \version 7.0.1 "Blackbird"
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
 

#include "../../include/toolboxes/C1DInterpolation.hpp"

/*--- Function for evaluating the value for akima spline ---*/
su2double CAkimaInterpolation::EvaluateSpline(su2double Point_Interp){
    Point_Match = true;
    int i =0, j=this->n-1;

    while (j-i>1){ 
        int m=(i+j) / 2 ; 
            if (Point_Interp>this->x[m]) {i=m}; 
            else {j=m}; 
    }

    su2double h=Point_Interp-this->x[i] ;
    return this->y[i]+h*(this->b[i]+h*(this->c[i]+h*this->d[i])) ;
}

/*--- Function for evaluating the value for linear 'spline' ---*/
su2double CLinearInterpolation::EvaluateSpline(su2double Point_Interp){
    int size = this->x.size();

    for (int i=0;i<size-1;i++){
        if(Point_Interp>=this->x[i] && Point_Interp<=this->x[i+1]){
            this->Point_Match = true;
            return (Point_Interp-this->x[i])*dydx[i]+this->y[i];}
    }
    return 0;
}

/*--- Function for setting the cofficients for linear 'spline' ---*/
void CLinearInterpolation::SetSpline(vector<su2double> &X, vector<su2double> &Data){
    int n = X.size();
    vector<su2double> h (n-1);
    this->x.resize(n);
    this->y.resize(n);
    this->dydx.resize(n-1);

    for (int i=0; i<n ; i++){
        this->x[i] = X[i];
        this->y[i] = Data[i];
    }

    for (int i=0; i<n-1;i ++){h[i]=X[i+1]-X[i];}
    for (int i=0; i<n-1;i++) {this->dydx[i]=(Data[i+1]-Data[i]) / h [i];}
}


/*--- Function for setting the cofficients for akima spline ---*/
void CAkimaInterpolation::SetSpline (vector<su2double> &X,vector<su2double> &Data){

    int n = X.size();
    vector<su2double> h (n-1);
    vector<su2double> p (n-1);

    /*---calculating finite differences (h) and gradients (p) ---*/
    for (int i=0; i<n-1;i ++){h[i]=X[i+1]-X[i];}
    for (int i=0; i<n-1;i++) {p[i]=(Data[i+1]-Data[i])/h[i];} 

    /*---b,c,d are the akima spline's cofficient for the cubic equation---*/
    this->x.resize(n);
    this->y.resize(n);
    this->b.resize(n);
    this->c.resize(n-1);
    this->d.resize(n-1);
    this->n=n; 

    for (int i=0; i<n ; i++){
        this->x[i]=X[i]; 
        this->y[i]=Data[i];
    }

    this->b[0] = p[0] ; 
    this->b[1] =(p[0]+ p[1])/2 ;
    this->b[n-1]=p[n-2]; 
    this->b [n-2]=(p[n-2]+p[n-3])/2;

    for (int i =2; i<n-2; i ++){
        su2double w1=fabs(p[i+1]-p[i]) , w2=fabs(p[i-1]-p[i-2]);

        if (w1+w2==0) {
            this->b[i] = (p[i-1]+p[i])/2 ;
        }
        else{
            this->b [i] = (w1*p[i-1]+w2*p[i])/(w1+w2);
        }
    }

    for (int i =0; i<n-1; i ++){
        this->c [i]=(3*p[i]-2*this->b[i]-this->b[i+1])/h[i] ;
        this->d [i]=(this->b[i+1]+this->b[i]-2*p[i])/h[i]/h[i];
    }
}

/*--- Correcting for Interpolation Type ---*/
vector<su2double> C1DInterpolation::CorrectedInletValues(vector<su2double> &Inlet_Interpolated , 
                                                        su2double Theta ,
                                                        unsigned short nDim, 
                                                        su2double *Coord, 
                                                        unsigned short nVar_Turb,
                                                        CConfig *config){
    su2double size_columns=Inlet_Interpolated.size()+nDim;
    vector<su2double> Inlet_Values (size_columns);
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
    switch(config->GetKindInletInterpolationType()){
    case(VR_VTHETA):
        unit_r = Inlet_Interpolated[nDim+2];
        unit_Theta = Inlet_Interpolated[nDim+3];
    break;
    case(ALPHA_PHI):
        Alpha = Inlet_Interpolated[nDim+2];
        Phi = Inlet_Interpolated[nDim+3];
        unit_m = sqrt(1/(1+pow(tan(Alpha),2)));
        unit_Theta = tan(Alpha)*unit_m;
        unit_r=unit_m*sin(Phi);
    break;
    }

    Inlet_Values[nDim+2] = unit_r*cos(Theta) - unit_Theta*sin(Theta); //for ix
    Inlet_Values[nDim+3] = unit_r*sin(Theta) + unit_Theta*cos(Theta); //for iy
    Inlet_Values[nDim+4] = sqrt(1-pow(unit_r,2)- pow(unit_Theta,2));  //for iz

    return Inlet_Values;

}

/*--- For printing interpolated data to a file ---*/
void C1DInterpolation::PrintInletInterpolatedData(vector<su2double>& Inlet_Values, string Marker, unsigned long nVertex, unsigned short nDim){
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

void C1DInterpolation::SetDataFromInletColumns(vector<su2double>& Inlet_Data, unsigned short nColumns, unsigned long nRows, unsigned short iCol, vector<su2double> &X, vector<su2double> &Data){
    int n = Inlet_Data.size();
    X.resize(n);
    Data.resize(n);
    
    unsigned long index;
    for (unsigned long iRow = 0; iRow < nRows; iRow++){
    index = iRow*nColumns;
        X[iRow]=Inlet_Data[index];
        Data[iRow]=Inlet_Data[index+iCol];
    }

    SetSpline(X,Data);
}