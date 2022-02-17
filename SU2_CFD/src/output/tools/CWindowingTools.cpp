/*!
 * \file signal_processing_toolbox.cpp
 * \brief Signal processing tools
 * \author S. Schotthöfer
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/output/tools/CWindowingTools.hpp"
#include <cassert>

//WindowingTools
su2double CWindowingTools::GetWndWeight(WINDOW_FUNCTION windowId, unsigned long curTimeIter, unsigned long endTimeIter) {
  switch (windowId) {
    case WINDOW_FUNCTION::HANN:        return HannWindow(curTimeIter, endTimeIter);
    case WINDOW_FUNCTION::HANN_SQUARE: return HannSquaredWindow(curTimeIter, endTimeIter);
    case WINDOW_FUNCTION::BUMP:        return BumpWindow(curTimeIter, endTimeIter);
    case WINDOW_FUNCTION::SQUARE:      return 1.0;
  }
  return 1.0;
}

su2double CWindowingTools::HannWindow(unsigned long curTimeIter, unsigned long endTimeIter) {
  su2double currTimeDouble = static_cast<su2double>(curTimeIter);
  if(endTimeIter==0) return 0; //Catch div by zero error, if window length is zero
  su2double endTimeDouble = static_cast<su2double>(endTimeIter);
  su2double tau = currTimeDouble/endTimeDouble;
  return 1.0-cos(2*PI_NUMBER*tau);
}

su2double CWindowingTools::HannSquaredWindow(unsigned long curTimeIter, unsigned long endTimeIter) {
  su2double currTimeDouble = static_cast<su2double>(curTimeIter);
  if(endTimeIter==0) return 0; //Catch div by zero error, if window length is zero
  su2double endTimeDouble = static_cast<su2double>(endTimeIter);
  su2double tau = currTimeDouble/endTimeDouble;
  return 2.0/3.0*(1-cos(2*PI_NUMBER*tau))*(1-cos(2*PI_NUMBER*tau));
}

su2double CWindowingTools::BumpWindow(unsigned long curTimeIter, unsigned long endTimeIter) {
  if(curTimeIter==0) return 0;
  if(curTimeIter==endTimeIter) return 0;
  su2double currTimeDouble = static_cast<su2double>(curTimeIter);
  su2double endTimeDouble = static_cast<su2double>(endTimeIter);
  su2double tau = currTimeDouble/endTimeDouble;
  return 1.0/0.00702986*(exp(-1/(tau-tau*tau)));
  /* 0.00702986 equals the integral of exp(-1/(tau-tau*tau)) from 0 to 1,
   * and it acts as a normalization constant */
}

void CWindowedAverage::addValue(su2double valIn, unsigned long curTimeIter,unsigned long startIter){
  su2double totalSum = 0;
  unsigned long windowWidth = 0;
  if (curTimeIter < startIter) return;  // Averaging not yet started.
  if (curTimeIter != lastTimeIter) {           // Handle new timestep
    if (curTimeIter > startIter) {
      windowWidth = curTimeIter - startIter;    // Calculate total width of window for this iteration
      cachedSum = updateCachedSum(windowWidth-1);  // Save weighted sum from last time step for later use
    } 
    lastTimeIter = curTimeIter;                // New time iteration step, update iteration number.
    // Add new sample
    if (windowingFunctionId != WINDOW_FUNCTION::SQUARE) {
      values.push_back(valIn);  // Add new sample to list for non-trivial windows
    }
  }
  else {  // We are within the same timestep. Update the last sample
    values.back() = valIn;
  }
  // Update the windowed-average from the weighted sum of previous samples and the latest sample
  totalSum = cachedSum + valIn*GetWndWeight(windowingFunctionId, windowWidth, windowWidth);
  val = totalSum / static_cast<su2double>(windowWidth);
}

su2double CWindowedAverage::updateCachedSum(unsigned long windowWidth) const { 
    su2double weightedSum = 0.;
    // Handle square window
    if (WINDOW_FUNCTION::SQUARE) return val * static_cast<su2double>(windowWidth);
    // Handle non-trivial windows
    if (values.size() == 0) return 0.;  // Handle first timestep
      for (unsigned long curTimeIter = 0; curTimeIter < values.size(); curTimeIter++) {
      weightedSum += values[curTimeIter] * GetWndWeight(windowingFunctionId, curTimeIter, values.size() - 1);
      }
      return weightedSum;
}


