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

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "../config_structure.hpp"
#include "../../SU2_CFD/include/CMarkerProfileReaderFVM.hpp"

using namespace std;

class CInletInterpolation{

    protected:
    vector<passivedouble> InletInterpolatedData, Inlet_Data;
    vector<su2double> Inlet_Values;
    unsigned short nDim;
    unsigned long nColumns, nRow;


    
    su2double LinearInterpolation(su2double Interp_Radius, unsigned long index,unsigned long iRow); 
    su2double AkimaInterpolation(su2double Interp_Radius, unsigned long index, unsigned long iRow_Akima); 
    su2double Interpolate(su2double Interp_Radius, unsigned long index,unsigned long iRow, CConfig *config);

    void Driver(CMarkerProfileReaderFVM profileReader, CGeometry *geometry, CConfig *config, unsigned short jMarker,unsigned short iMarker);
    void SetInterpolatedData(unsigned long iVertex); 
    void PrintInterpolatedData(CMarkerProfileReaderFVM profileReader, CGeometry *geometry, unsigned short jMarker,unsigned short iMarker);
    void CorrectForInterpolationType(su2double Parameter1, su2double Parameter2, CConfig *config, su2double Theta);

    su2double Get_Ai_dash(unsigned long index, unsigned long iRow_Akima);
    su2double Get_Pi(unsigned long index, unsigned long iRow_Akima);
    su2double Get_Wi(unsigned long index, unsigned long iRow_Akima);

    bool Point_Match = false;

    public:
    
    CInletInterpolation(CGeometry *geometry, CConfig *config,CMarkerProfileReaderFVM profileReader, string profile_filename, unsigned short KIND_MARKER,unsigned short iMarker, unsigned short jMarker, unsigned short nDim);
    ~CInletInterpolation(void);

    inline vector<passivedouble> GetInterpolatedProfile()
    {return InletInterpolatedData;}

};
