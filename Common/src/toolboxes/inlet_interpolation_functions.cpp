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
 

#include "../../include/toolboxes/inlet_interpolation_functions.hpp"

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
    bi = Get_Ai_dash(Inlet_Data, iRow, index,nColumns,nRow);
    ci = (3*Get_Pi(Inlet_Data, iRow, index,nColumns)-2*bi-Get_Ai_dash(Inlet_Data, iRow+1, index, nColumns, nRow))/dxi;
    di = (bi + Get_Ai_dash(Inlet_Data, iRow+1, index,nColumns, nRow) - 2*Get_Pi(Inlet_Data, iRow,index,nColumns))/pow(dxi,2);
    delta = Interp_Radius - Inlet_Data[nColumns*iRow];    

    interpolated_value = ai+bi*delta+ci*pow(delta,2)+di*pow(delta,3);

    return interpolated_value;

}

su2double CSolver::Get_Ai_dash(vector<su2double> Inlet_Data, unsigned long iRow, unsigned short index, unsigned short nColumns, unsigned long nRow){
    //for Boundary conditions (two first and two last points require special definition, can vary for different codes)
    if(iRow == 0) {return Get_Pi(Inlet_Data, iRow,index,nColumns);}
    else if (iRow == 1) {return (Get_Pi(Inlet_Data, iRow-1,index,nColumns)+Get_Pi(Inlet_Data, iRow,index,nColumns))/2;}
    else if (iRow == nRow-2) {return (Get_Pi(Inlet_Data,nRow-2,index,nColumns)+Get_Pi(Inlet_Data, nRow-3,index, nColumns))/2;}
    else if (iRow == nRow-1) {return Get_Pi(Inlet_Data, nRow-2,index,nColumns);}
    else{
    if((Get_Wi(Inlet_Data, iRow+1,index, nColumns)+Get_Wi(Inlet_Data, iRow-1,index,nColumns)) != 0)
      return (Get_Wi(Inlet_Data, iRow+1,index,nColumns)*Get_Pi(Inlet_Data, iRow-1,index,nColumns) + Get_Wi(Inlet_Data, iRow-1,index,nColumns)*Get_Pi(Inlet_Data, iRow,index,nColumns))/(Get_Wi(Inlet_Data, iRow+1,index,nColumns) + Get_Wi(Inlet_Data, iRow-1,index,nColumns));
    else
      return ((Get_Pi(Inlet_Data, iRow-1,index,nColumns)-Get_Pi(Inlet_Data, iRow,index,nColumns))/2);
    }
}


su2double CSolver::Get_Pi(vector<su2double> Inlet_Data, unsigned long iRow, unsigned short index, unsigned short nColumns) {
return (Inlet_Data[nColumns*(iRow+1)+index] - Inlet_Data[nColumns*iRow+index])/(Inlet_Data[nColumns*(iRow+1)] - Inlet_Data[nColumns*iRow]);}

su2double CSolver::Get_Wi(vector<su2double> Inlet_Data, unsigned long iRow,unsigned short index,unsigned short nColumns) {return fabs(Get_Pi(Inlet_Data,iRow,index,nColumns) - Get_Pi(Inlet_Data,iRow-1,index,nColumns));}
