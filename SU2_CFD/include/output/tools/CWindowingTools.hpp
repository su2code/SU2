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

#include <vector>
#include "../../../Common/include/option_structure.hpp"

class CWindowingTools{
public:

    /*! \brief Returns the value of a windowing function given by fctIdx at CurTimeIdx with given endTimeIdx (i.e. endTimeIdx=nTimeIter, if one starts  windowing at time t =0.) */
    su2double GetWndWeight(WINDOW_FUNCTION windowId, unsigned long CurTimeIdx, unsigned long endTimeIdx);

protected:
// Long time windows
    su2double HannWindow(unsigned long i, unsigned long endTimeIdx);

    su2double HannSquaredWindow(unsigned long i, unsigned long endTimeIdx);

    su2double BumpWindow(unsigned long i, unsigned long endTimeIdx);
};

class CWindowedAverage:CWindowingTools{
private:
su2double val;
unsigned long count;
std::vector<su2double> values;


public:
CWindowedAverage();

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
