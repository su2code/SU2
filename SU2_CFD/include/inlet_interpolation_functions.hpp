#pragma once

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "../../Common/include/config_structure.hpp"
//#include "CMarkerProfileReaderFVM.hpp"

using namespace std;

su2double CSolver::ONEDLINEAR_SPANWISE(vector<su2double> Inlet_Data, su2double Interp_Radius, unsigned long iRow,unsigned short index,unsigned long nColumns)
{
    su2double slope, interpolated_value;

    slope = (Inlet_Data[index+nColumns*(iRow+1)]-Inlet_Data[index+nColumns*iRow])/(Inlet_Data[nColumns*(iRow+1)]-Inlet_Data[nColumns*iRow]);

    interpolated_value = slope*(Interp_Radius - Inlet_Data[nColumns*iRow]) + Inlet_Data[index+nColumns*iRow];

    return interpolated_value;
}

su2double CSolver::ONEDAKIMA_SPANWISE(vector<su2double> Inlet_Data, su2double Interp_Radius, unsigned long iRow,unsigned short index,unsigned short nColumns, unsigned long nRow)
{
    su2double dxi, ai, bi, ci ,di ,delta, interpolated_value;
    /*--- Finding the cofficients of the third order polynomial for Akima Interpolation ---*/
    dxi = Inlet_Data[nColumns*(iRow+1)] - Inlet_Data[nColumns*iRow];
    ai = Inlet_Data[nColumns*iRow + index];
    bi = Get_Ai_dash(Inlet_Data, iRow, index, nRow);
    ci = (3*Get_Pi(Inlet_Data, iRow, index)-2*bi-Get_Ai_dash(Inlet_Data, iRow+1, index, nRow))/dxi;
    di = (bi + Get_Ai_dash(Inlet_Data, iRow+1, index, nRow) - 2*Get_Pi(Inlet_Data, iRow,index))/pow(dxi,2);
    delta = Interp_Radius - Inlet_Data[nColumns*iRow];    

    interpolated_value = ai+bi*delta+ci*pow(delta,2)+di*pow(delta,3);

    return interpolated_value;

}

su2double CSolver::Get_Ai_dash(vector<su2double> Inlet_Data, unsigned long iRow, unsigned short index, unsigned short nColumns, unsigned long nRow){
    //for Boundary conditions (two first and two last points require special definition, can vary for different codes)
    if(iRow == 0) {return Get_Pi(Inlet_Data, iRow,index);}
    else if (iRow == 1) {return (Get_Pi(Inlet_Data, iRow-1,index)+Get_Pi(Inlet_Data, iRow,index))/2;}
    else if (iRow == nRow-2) {return (Get_Pi(Inlet_Data,nRow-2,index)+Get_Pi(Inlet_Data, nRow-3,index))/2;}
    else if (iRow == nRow-1) {return Get_Pi(Inlet_Data, nRow-2,index);}
    else{
    if((Get_Wi(Inlet_Data, iRow+1,index)+Get_Wi(Inlet_Data, iRow-1,index)) != 0)
      return (Get_Wi(Inlet_Data, iRow+1,index)*Get_Pi(Inlet_Data, iRow-1,index) + Get_Wi(Inlet_Data, iRow-1,index)*Get_Pi(Inlet_Data, iRow,index))/(Get_Wi(Inlet_Data, iRow+1,index) + Get_Wi(Inlet_Data, iRow-1,index));
    else
      return ((Get_Pi(Inlet_Data, iRow-1,index)-Get_Pi(Inlet_Data, iRow,index))/2);
    }
}


su2double CSolver::Get_Pi(vector<su2double> Inlet_Data, unsigned long iRow, unsigned short index, unsigned short nColumns) {
return (Inlet_Data[nColumns*(iRow+1)+index] - Inlet_Data[nColumns*iRow+index])/(Inlet_Data[nColumns*(iRow+1)] - Inlet_Data[nColumns*iRow]);}

su2double CSolver::Get_Wi(vector<su2double> Inlet_Data, unsigned long iRow,unsigned short index) {return fabs(Get_Pi(Inlet_Data,iRow,index) - Get_Pi(Inlet_Data,iRow-1,index));}
