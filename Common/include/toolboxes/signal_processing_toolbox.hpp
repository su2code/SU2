#pragma once

#include "../datatype_structure.hpp"
#include <vector>

namespace Signal_Processing {
  
  su2double Average(std::vector<su2double> &data);
  
  class RunningAverage
  {
  private:
    su2double val;
    unsigned long count;
    const su2double PI_NUMBER = 4.0 * atan(1.0);	/*!< \brief Pi number. */

  
  public:
    RunningAverage(){
      this->Reset();
    }
 
    su2double Update(su2double valIn){
      su2double scaling = 1. / (su2double)(count + 1);
      val = valIn * scaling + val * (1. - scaling);
      count++;
      return val;
    }
  
    su2double Get(){
      return val;
    }
  
    unsigned long Count(){
      return count;
    }
  
    void Reset(){
      val = 0.;
      count = 0;
    }

    /*! \brief Returns the value of a windowing function given by fctIdx at CurTimeIdx with given endTimeIdx (i.e. endTimeIdx=nTimeIter, if one starts  windowing at time t =0.) */
    su2double GetWndWeight(int fctIdx, unsigned long CurTimeIdx, unsigned long endTimeIdx){
        switch (fctIdx){
          case 1: return HannWindow(CurTimeIdx, endTimeIdx);
          case 2: return HannSquaredWindow(CurTimeIdx, endTimeIdx);
          case 3: return BumpWindow(CurTimeIdx, endTimeIdx);
          default:return 1.0;
        }
    }
    private:
    su2double HannWindow(unsigned long i, unsigned long endTimeIdx){
      if(i==0) return 0; //catch div by zero error
      su2double currTime = static_cast<su2double>(i);
      su2double endTime = static_cast<su2double>(endTimeIdx);
      su2double tau = currTime/endTime;
      return 1.0-cos(2*PI_NUMBER*tau);
    }

    su2double HannSquaredWindow(unsigned long i, unsigned long endTimeIdx){
      if(i==0) return 0; //catch div by zero error
      su2double currTime = static_cast<su2double>(i);
      su2double endTime = static_cast<su2double>(endTimeIdx);
      su2double tau = currTime/endTime;
      return 2.0/3.0*(1-cos(2*PI_NUMBER*tau))*(1-cos(2*PI_NUMBER*tau));
    }

    su2double BumpWindow(unsigned long i, unsigned long endTimeIdx){
      if(i==0) return 0;
      if(i==endTimeIdx) return 0;
      su2double currTime = static_cast<su2double>(i);
      su2double endTime = static_cast<su2double>(endTimeIdx);
      su2double tau = currTime/endTime;
      return 1.0/0.00702986*(exp(-1/(tau-tau*tau)));
    }
  };
};
