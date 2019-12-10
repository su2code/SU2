/*!
 * \file signal_processing_toolbox.hpp
 * \brief Header file for the signal processing toolbox.
 * \author T. Albring
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

#include "../datatype_structure.hpp"
#include <vector>


class WindowingTools{
public:
    WindowingTools(){}
    ~WindowingTools(){}

    /*! \brief Returns the value of a windowing function given by fctIdx at CurTimeIdx with given endTimeIdx (i.e. endTimeIdx=nTimeIter, if one starts  windowing at time t =0.) */
    su2double GetWndWeight(int fctIdx, unsigned long CurTimeIdx, unsigned long endTimeIdx);

protected:
    const su2double PI_NUMBER = 4.0 * atan(1.0);	/*!< \brief Pi number. */

// Long time windows
    su2double HannWindow(unsigned long i, unsigned long endTimeIdx);

    su2double HannSquaredWindow(unsigned long i, unsigned long endTimeIdx);

    su2double BumpWindow(unsigned long i, unsigned long endTimeIdx);

// Short time windows
    /*su2double TriangleWindow(unsigned long i, unsigned long endTimeIdx){
        su2double currTime = static_cast<su2double>(i);
        su2double endTime = static_cast<su2double>(endTimeIdx);
        su2double tau = currTime/endTime;
        if(tau < 0.5) return 4*tau;
        else          return 4*(1-tau);
    }

    su2double QuadWindow(unsigned long i, unsigned long endTimeIdx){
        su2double currTime = static_cast<su2double>(i);
        su2double endTime = static_cast<su2double>(endTimeIdx);
        su2double tau = currTime/endTime;
        if(tau < 1.0/3.0) return 13.5*tau*tau;
        if(tau < 2.0/3.0) return 27.0*tau*(1-tau)-4.5;
        else              return 13.5*(1-tau)*(1-tau);
    }
    */
};

class WindowedAverage:WindowingTools{
private:
su2double val;
unsigned long count;
std::vector<su2double> values;


public:
WindowedAverage();

~WindowedAverage(){}

su2double Update(su2double valIn);

su2double Get();

unsigned long Count();

void Reset();

void addValue(su2double valIn, unsigned long currentIter,unsigned long startIter = 0);

su2double WindowedUpdate(int fctIdx); //Computes a windowed time average (integral)

private:
//Using Midpoint rule for consistentcy with adjoint solver
su2double NoWindowing(); // == Square Window

su2double HannWindowing();

su2double HannSquaredWindowing();

su2double BumpWindowing();
};
