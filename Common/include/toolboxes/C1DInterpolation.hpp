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
    bool Point_Match = true;
public:
virtual void SetSpline(vector<su2double> &x, vector<su2double> &y){}
virtual su2double EvaluateSpline(su2double Point_Interp){}
virtual bool GetPointMatch(){return Point_Match;}
/*
void Interpolation_Set(vector<su2double> &x, vector<su2double> &y, CConfig *config){
switch(config->GetKindInletInterpolationFunction()){
    case (ONED_AKIMASPLINE_SPANWISE):
        return this->SetAkimaSpline(x,y);
    case(ONED_LINEAR_SPANWISE):
        return this->SetLinearSpline(x,y);
}

}

su2double Interpolation_Evaluate(su2double Point_Interp, string Interpolation_Function){
        if(Interpolation_Type == "Akima")
        return this->EvalAkimaSpline(Point_Interp);
    else if(Interpolation_Type == "Linear")
        return this->EvalLinearSpline(Point_Interp);
}
*/
vector<su2double> CorrectedInletValues(vector<su2double> &Inlet_Interpolated, 
                                    su2double Theta ,
                                    unsigned short nDim, 
                                    su2double *Coord, 
                                    unsigned short nVar_Turb, 
                                    CConfig *config);

void PrintInletInterpolatedData(vector<su2double> Inlet_Values, string Marker, unsigned long nVertex, unsigned short nDim);
};


class CAkimaInterpolation: public C1DInterpolation{ 
protected:
    vector<su2double> x,y,b,c,d;
    int n;
public:
/*   CAkimaInterpolation(vector<su2double> &x, vector<su2double> &y, unsigned short nColumns){
        SetAkimaSpline(x, y);


}*/
    void SetSpline(vector<su2double> &x, vector<su2double> &y);
    su2double EvaluateSpline(su2double Point_Interp);
    bool GetPointMatch(){return Point_Match;}
};

class CLinearInterpolation: public C1DInterpolation{
    protected:
    vector<su2double> x,y,dydx;
    public:
    void SetSpline(vector<su2double> &x, vector<su2double> &y);
    su2double EvaluateSpline(su2double Point_Interp);
    //bool Point_Match = false;
    bool GetPointMatch(){return Point_Match;}
};
