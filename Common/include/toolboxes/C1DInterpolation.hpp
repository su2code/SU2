/*!
 * \file inlet_interpolation_functions.hpp
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

#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include<fstream>
#include "../datatype_structure.hpp"
#include "../CConfig.hpp"

using namespace std;

class C1DInterpolation{
protected:
    bool Point_Match = false;
    vector<su2double> Data;
    vector<su2double> X;
public:
virtual void SetSpline(vector<su2double> &x, vector<su2double> &y){}
virtual su2double EvaluateSpline(su2double Point_Interp){return 0;}
virtual bool GetPointMatch(){return Point_Match;}

vector<su2double> CorrectedInletValues(vector<su2double> &Inlet_Interpolated, 
                                    su2double Theta ,
                                    unsigned short nDim, 
                                    su2double *Coord, 
                                    unsigned short nVar_Turb, 
                                    CConfig *config);

void PrintInletInterpolatedData(vector<su2double> Inlet_Values, string Marker, unsigned long nVertex, unsigned short nDim);

void SetDataFromInletColumns(vector<su2double>& Inlet_Data, unsigned short nColumns, unsigned long nRows, unsigned short iCol);
};


class CAkimaInterpolation: public C1DInterpolation{ 
protected:
    vector<su2double> x,y,b,c,d;
    int n;
public:
    CAkimaInterpolation(vector<su2double>& Inlet_Data, unsigned short nColumns, unsigned long nRows, signed short iCol){
    
    SetDataFromInletColumns(Inlet_Data, nColumns, nRows, iCol);
    SetSpline(this->X, this->Data);
    }

    void SetSpline(vector<su2double> &x, vector<su2double> &y);
    su2double EvaluateSpline(su2double Point_Interp);
    bool GetPointMatch(){return Point_Match;}
};

class CLinearInterpolation: public C1DInterpolation{
    protected:
    vector<su2double> x,y,dydx;
    public:
    CLinearInterpolation(vector<su2double>& Inlet_Data, unsigned short nColumns, unsigned long nRows, signed short iCol){
    
    SetDataFromInletColumns(Inlet_Data, nColumns, nRows, iCol);
    SetSpline(this->X, this->Data);
    }

    void SetSpline(vector<su2double> &x, vector<su2double> &y);
    su2double EvaluateSpline(su2double Point_Interp);
    //bool Point_Match = false;
    bool GetPointMatch(){return Point_Match;}
};
