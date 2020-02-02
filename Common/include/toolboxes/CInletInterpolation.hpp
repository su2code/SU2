/*!
 * \file inlet_interpolation_functions.hpp
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

#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include<fstream>
#include "../datatype_structure.hpp"

using namespace std;

class CAkimaSpline{ 
protected:
vector<su2double> x,y,b,c,d;
int n;
public:
void akima_spline_set(CAkimaSpline *s, vector<su2double> &x, vector<su2double> &y);
su2double akima_spline_eval(CAkimaSpline *s, su2double Point_Interp);
};

class CLinearSpline{
    protected:
    vector<su2double> x,y,dydx;
    public:
    void Linear_set(CLinearSpline *s,vector<su2double> &x, vector<su2double> &y);
    su2double Linear_eval(CLinearSpline *s, su2double Point_Interp);
    bool Point_Match = false;
};

class CInletInterpolation:public CAkimaSpline,public CLinearSpline{
string Interpolation_Type, Interpolation_Function;
public:
void Interpolation_Set(CInletInterpolation *s, vector<su2double> &x, vector<su2double> &y, string Interpolation_Function){
    if(Interpolation_Type == "Akima")
        return s->akima_spline_set(s,x,y);
    else if(Interpolation_Type == "Linear")
        return s->Linear_set(s,x,y);
}

su2double Interpolation_Evaluate(CInletInterpolation *s, su2double Point_Interp, string Interpolation_Function){
        if(Interpolation_Type == "Akima")
        return s->akima_spline_eval(s,Point_Interp);
    else if(Interpolation_Type == "Linear")
        return s->Linear_eval(s,Point_Interp);
}

vector<su2double>& CorrectedInletValues(vector<su2double> &Inlet_Interpolated, 
                                    su2double Theta ,
                                    unsigned short nDim, 
                                    su2double *Coord, 
                                    unsigned short nVar_Turb, 
                                    string Interpolation_Type);

void PrintInletInterpolatedData(vector<su2double> Inlet_Values, string Marker, unsigned long nVertex, unsigned short nDim);

void free_memory (CInletInterpolation *s){
   delete[] s; s=NULL; }
};
