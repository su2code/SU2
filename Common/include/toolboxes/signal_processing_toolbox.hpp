#pragma once

#include "../datatype_structure.hpp"
#include <vector>

namespace Signal_Processing {
  
  su2double Average(const std::vector<su2double> &data);
  
  class RunningAverage
  {
  private:
    su2double val;
    unsigned long count;
  
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
  
    su2double Get() const{
      return val;
    }
  
    unsigned long Count() const{
      return count;
    }
  
    void Reset(){
      val = 0.;
      count = 0;
    }
  };
};
