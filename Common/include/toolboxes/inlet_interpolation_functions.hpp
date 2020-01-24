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

    CMarkerProfileReaderFVM profileReader(CGeometry **geometry, CConfig *config,string profile_filename,unsigned short KIND_MARKER, unsigned short iMarker, unsigned short jMarker, unsigned short nDim);

    unsigned long nColumns, iPoint, iVertex, iRow, index, iRow_Akima, nRow;
    unsigned short nCol_InletFile;
    unsigned short iMarker, jMarker;
    CGeometry **geometry;
    CConfig *config;

    su2double Interp_Radius, Theta, *Coord;

    su2double slope, interpolated_value, Parameter1, Parameter2, unit_r, unit_Theta, unit_m, Alpha, Phi;
    su2double dxi, ai, bi, ci ,di ,delta;
    
    void LinearInterpolation(), AkimaInterpolation(unsigned long iRow_Akima), CorrectForInterpolationType(), SetInterpolatedData(), PrintInterpolatedData(CMarkerProfileReaderFVM profileReader);
    void Interpolate();
    void GetClosestPointFromFile();
    void SetVertex(CMarkerProfileReaderFVM profileReader);

    su2double Get_Ai_dash(unsigned long iRow_Akima), Get_Pi(unsigned long iRow_Akima), Get_Wi(unsigned long iRow_Akima);

    bool Point_Match = false;

    public:
    
    CInletInterpolation(CGeometry **geometry, CConfig *config, string profile_filename, unsigned short KIND_MARKER,unsigned short iMarker, unsigned short jMarker, unsigned short nDim);
    ~CInletInterpolation(void);

    inline vector<passivedouble> GetInterpolatedProfile()
    {return InletInterpolatedData;}

    inline unsigned long GetNumberofColumns(){
        return nColumns+nDim;}

    inline unsigned long GetNumberofVertexes(){
        return geometry[MESH_0]->nVertex[iMarker];
    }

};
