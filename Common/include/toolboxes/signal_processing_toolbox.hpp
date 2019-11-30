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
