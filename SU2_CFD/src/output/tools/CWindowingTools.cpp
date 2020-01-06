/*!
 * \file signal_processing_toolbox.cpp
 * \brief Signal processing tools
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

#include "../../../include/output/tools/CWindowingTools.hpp"

//WindowingTools
/*! \brief Returns the value of a windowing function given by fctIdx at CurTimeIdx with given endTimeIdx (i.e. endTimeIdx=nTimeIter, if one starts  windowing at time t =0.) */
su2double CWindowingTools::GetWndWeight(WINDOW_FUNCTION windowId, unsigned long CurTimeIdx, unsigned long endTimeIdx) const{
  switch (windowId) {
    case HANN:        return HannWindow(CurTimeIdx, endTimeIdx);
    case HANN_SQUARE: return HannSquaredWindow(CurTimeIdx, endTimeIdx);
    case BUMP:        return BumpWindow(CurTimeIdx, endTimeIdx);
    default:return 1.0;
  }
}
su2double CWindowingTools::HannWindow(unsigned long i, unsigned long endTimeIdx) const{
  su2double currTime = static_cast<su2double>(i);
  if(endTimeIdx==0) return 0; //Catch div by zero error, if window length is zero
  su2double endTime = static_cast<su2double>(endTimeIdx);
  su2double tau = currTime/endTime;
  return 1.0-cos(2*PI_NUMBER*tau);
}

su2double CWindowingTools::HannSquaredWindow(unsigned long i, unsigned long endTimeIdx) const{
  su2double currTime = static_cast<su2double>(i);
  if(endTimeIdx==0) return 0; //Catch div by zero error, if window length is zero
  su2double endTime = static_cast<su2double>(endTimeIdx);
  su2double tau = currTime/endTime;
  return 2.0/3.0*(1-cos(2*PI_NUMBER*tau))*(1-cos(2*PI_NUMBER*tau));
}

su2double CWindowingTools::BumpWindow(unsigned long i, unsigned long endTimeIdx) const{
  if(i==0) return 0;
  if(i==endTimeIdx) return 0;
  su2double currTime = static_cast<su2double>(i);
  su2double endTime = static_cast<su2double>(endTimeIdx);
  su2double tau = currTime/endTime;
  return 1.0/0.00702986*(exp(-1/(tau-tau*tau)));
  /* 0.00702986 equals the integral of exp(-1/(tau-tau*tau)) from 0 to 1,
   * and it acts as a normalization constant */
}


//WindowedAverage
CWindowedAverage::CWindowedAverage(){
  this->Reset();
}

su2double CWindowedAverage::Update(su2double valIn){
  su2double scaling = 1. / (su2double)(count + 1);
  val = valIn * scaling + val * (1. - scaling);
  count++;
  return val;
}

su2double CWindowedAverage::Get() const{
  return val;
}

unsigned long CWindowedAverage::Count() const{
  return count;
}

void CWindowedAverage::Reset(){
  val = 0.;
  count = 0;
}


void CWindowedAverage::addValue(su2double valIn, unsigned long currentIter,unsigned long startIter){
  if(currentIter >= startIter)values.push_back(valIn);
}

su2double CWindowedAverage::WindowedUpdate(int fctIdx){ //Computes a windowed time average (integral)
  if(values.size()>0){
    switch (fctIdx){
      case HANN:        val= HannWindowing();         return val;
      case HANN_SQUARE: val= HannSquaredWindowing();  return val;
      case BUMP:        val= BumpWindowing();         return val;
      default:          val= NoWindowing();           return val;
    }
  }
  return 0.0;
}

su2double CWindowedAverage::NoWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
    wnd_timeAvg+=values[i];
  }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}

su2double CWindowedAverage::HannWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
    wnd_timeAvg+=values[i]*HannWindow(i,values.size()-1);
  }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}

su2double CWindowedAverage::HannSquaredWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
    wnd_timeAvg+=values[i]*HannSquaredWindow(i,values.size()-1);
  }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}

su2double CWindowedAverage::BumpWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
    wnd_timeAvg+=values[i]*BumpWindow(i,values.size()-1);
  }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}


