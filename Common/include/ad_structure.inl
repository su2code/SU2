/*!
 * \file ad_structure.inl
 * \brief Main routines for the algorithmic differentiation (AD) structure.
 * \author T. Albring
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

namespace AD{
#if defined CODI_REVERSE_TYPE

  typedef codi::DataStore CheckpointHandler;

  typedef su2double::TapeType Tape;

  typedef codi::ExternalFunctionHelper<su2double> ExtFuncHelper;

  extern ExtFuncHelper* FuncHelper;

  /*--- Stores the indices of the input variables (they might be overwritten) ---*/

  extern std::vector<su2double::GradientData> inputValues;

  /*--- Current position inside the adjoint vector ---*/

  extern int adjointVectorPosition;

  /*--- Reference to the tape ---*/

  extern su2double::TapeType& globalTape;

  extern bool Status;

  extern bool PreaccActive;

  extern bool PreaccEnabled;

  extern su2double::TapeType::Position StartPosition, EndPosition;

  extern std::vector<su2double::TapeType::Position> TapePositions;

  extern std::vector<su2double::GradientData> localInputValues;

  extern std::vector<su2double*> localOutputValues;

  extern codi::PreaccumulationHelper<su2double> PreaccHelper;  

  inline void RegisterInput(su2double &data, bool push_index) {
    AD::globalTape.registerInput(data);
    if (push_index) {
      inputValues.push_back(data.getGradientData());
    }
  }

  inline void RegisterOutput(su2double& data) {AD::globalTape.registerOutput(data);}

  inline void ResetInput(su2double &data) {data.getGradientData() = su2double::GradientData();}

  inline void StartRecording() {AD::globalTape.setActive();}

  inline void StopRecording() {AD::globalTape.setPassive();}

  inline bool TapeActive() { return AD::globalTape.isActive(); }

  inline void PrintStatistics() {AD::globalTape.printStatistics();}

  inline void ClearAdjoints() {AD::globalTape.clearAdjoints(); }

  inline void ComputeAdjoint() {AD::globalTape.evaluate();
                               adjointVectorPosition = 0;}

  inline void ComputeAdjoint(unsigned short enter, unsigned short leave) {
    AD::globalTape.evaluate(TapePositions[enter], TapePositions[leave]);
    if (leave == 0) {
      adjointVectorPosition = 0;
    }
  }

  inline void Reset() {
    globalTape.reset();
    if (inputValues.size() != 0) {
      adjointVectorPosition = 0;
      inputValues.clear();
    }
    if (TapePositions.size() != 0) {
      TapePositions.clear();
    }    
  }

  inline void SetIndex(int &index, const su2double &data) {
    index = data.getGradientData();
  }

  inline void SetDerivative(int index, const double val) {
    AD::globalTape.setGradient(index, val);
  }

  inline double GetDerivative(int index) {
    return AD::globalTape.getGradient(index);
  }

  inline void SetPreaccIn(const su2double &data) {
    if (PreaccActive) {
      if (data.isActive()) {
        PreaccHelper.addInput(data);       
      }
    }
  }

  inline void SetPreaccIn(const su2double* data, const int size) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size; i++) {
        if (data[i].isActive()) {
          PreaccHelper.addInput(data[i]);
        }
      }
    }
  }

  inline void SetPreaccIn(const su2double* const *data, const int size_x, const int size_y) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size_x; i++) {
        for (unsigned short j = 0; j < size_y; j++) {
          if (data[i][j].isActive()) {
            PreaccHelper.addInput(data[i][j]);
          }
        }
      }
    }
  }

  inline void StartPreacc() {
    if (globalTape.isActive() && PreaccEnabled) {
      PreaccHelper.start();
      PreaccActive = true;
    }
  }

  inline void SetPreaccOut(su2double& data) {
    if (PreaccActive) {
      if (data.isActive()) {
        PreaccHelper.addOutput(data);
      }
    }
  }

  inline void SetPreaccOut(su2double* data, const int size) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size; i++) {
        if (data[i].isActive()) {
          PreaccHelper.addOutput(data[i]);
        }
      }
    }
  }

  inline void SetPreaccOut(su2double** data, const int size_x, const int size_y) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size_x; i++) {
        for (unsigned short j = 0; j < size_y; j++) {
          if (data[i][j].isActive()) {
            PreaccHelper.addOutput(data[i][j]);
          }
        }
      }
    }
  }

  inline void Push_TapePosition() {
    TapePositions.push_back(AD::globalTape.getPosition());
  }

  inline void EndPreacc(){
    if (PreaccActive) {
      PreaccHelper.finish(false);
    }
  }
  
  inline void StartExtFunc(bool storePrimalInput, bool storePrimalOutput){
    FuncHelper = new ExtFuncHelper(true);
    if (!storePrimalInput){
      FuncHelper->disableInputPrimalStore();
    }
    if (!storePrimalOutput){
      FuncHelper->disableOutputPrimalStore();
    }
  }
  
  inline void SetExtFuncIn(const su2double &data) {
    FuncHelper->addInput(data);       
  }

  inline void SetExtFuncIn(const su2double* data, const int size) {
    for (int i = 0; i < size; i++) {
      FuncHelper->addInput(data[i]);
    }

  }

  inline void SetExtFuncIn(const su2double* const *data, const int size_x, const int size_y) {
    for (int i = 0; i < size_x; i++) {
      for (int j = 0; j < size_y; j++) {
        FuncHelper->addInput(data[i][j]);
      }
    }
  }
  
  inline void SetExtFuncOut(su2double& data) {
    if (globalTape.isActive()) {
      FuncHelper->addOutput(data);
    }
  }

  inline void SetExtFuncOut(su2double* data, const int size) {
    for (int i = 0; i < size; i++) {
      if (globalTape.isActive()) {
        FuncHelper->addOutput(data[i]);
      }
    }
  }

  inline void SetExtFuncOut(su2double** data, const int size_x, const int size_y) {
    for (int i = 0; i < size_x; i++) {
      for (int j = 0; j < size_y; j++) {
        if (globalTape.isActive()) {
          FuncHelper->addOutput(data[i][j]);
        }
      }
    }
  }

  inline void delete_handler(void *handler) {
    CheckpointHandler *checkpoint = static_cast<CheckpointHandler*>(handler);
    checkpoint->clear();
  }
  
  inline void EndExtFunc(){delete FuncHelper;}
  
#else

  /*--- Default implementation if reverse mode is disabled ---*/

  inline void RegisterInput(su2double &data, bool push_index) {}

  inline void RegisterOutput(su2double& data) {}

  inline void StartRecording() {}

  inline void StopRecording() {}

  inline bool TapeActive() { return false; }

  inline void PrintStatistics() {}

  inline void ClearAdjoints() {}

  inline void ComputeAdjoint() {}

  inline void ComputeAdjoint(unsigned short enter, unsigned short leave) {}

  inline void SetIndex(int &index, const su2double &data) {}

  inline void SetDerivative(int index, const double val) {}

  inline double GetDerivative(int position) { return 0.0; }

  inline void Reset() {}

  inline void ResetInput(su2double &data) {}

  inline void SetPreaccIn(const su2double &data) {}

  inline void SetPreaccIn(const su2double* data, const int size) {}

  inline void SetPreaccIn(const su2double* const *data, const int size_x, const int size_y) {}

  inline void SetPreaccOut(su2double &data) {}

  inline void SetPreaccOut(su2double* data, const int size) {}

  inline void SetPreaccOut(su2double** data, const int size_x, const int size_y) {}

  inline void StartPreacc() {}

  inline void EndPreacc() {}

  inline void Push_TapePosition() {}

  inline void StartExtFunc(bool storePrimalInput, bool storePrimalOutput){}
  
  inline void SetExtFuncIn(const su2double &data) {}

  inline void SetExtFuncIn(const su2double* data, const int size) {}

  inline void SetExtFuncIn(const su2double* const *data, const int size_x, const int size_y) {}
  
  inline void SetExtFuncOut(su2double& data) {}

  inline void SetExtFuncOut(su2double* data, const int size) {}

  inline void SetExtFuncOut(su2double** data, const int size_x, const int size_y) {}
  
  inline void EndExtFunc(){}
#endif
}

/*--- If we compile under OSX we have to overload some of the operators for
 *   complex numbers to avoid the use of the standard operators
 *  (they use a lot of functions that are only defined for doubles) ---*/

#ifdef __APPLE__

namespace std{
  
  template<>
  inline su2double abs(const complex<su2double>& x){
    
    return sqrt(x.real()*x.real() + x.imag()*x.imag());
    
  }
  
  template<>
  inline complex<su2double> operator/(const complex<su2double>& x,
                                      const complex<su2double>& y){
    
    su2double d    = (y.real()*y.real() + y.imag()*y.imag());
    su2double real = (x.real()*y.real() + x.imag()*y.imag())/d;
    su2double imag = (x.imag()*y.real() - x.real()*y.imag())/d;
    
    return complex<su2double>(real, imag);
    
  }
  
  template<>
  inline complex<su2double> operator*(const complex<su2double>& x,
                                      const complex<su2double>& y){
    
    su2double real = (x.real()*y.real() - x.imag()*y.imag());
    su2double imag = (x.imag()*y.real() + x.real()*y.imag());
    
    return complex<su2double>(real, imag);
    
  }
}
#endif

