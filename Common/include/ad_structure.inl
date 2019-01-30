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

  extern std::vector<su2double::GradientData> localInputValues;

  extern std::vector<su2double*> localOutputValues;

  inline void RegisterInput(su2double &data) {AD::globalTape.registerInput(data);
                                             inputValues.push_back(data.getGradientData());}

  inline void RegisterOutput(su2double& data) {AD::globalTape.registerOutput(data);}

  inline void ResetInput(su2double &data) {data.getGradientData() = su2double::GradientData();}

  inline void StartRecording() {AD::globalTape.setActive();}

  inline void StopRecording() {AD::globalTape.setPassive();}

  inline void ClearAdjoints() {AD::globalTape.clearAdjoints(); }

  inline void ComputeAdjoint() {AD::globalTape.evaluate();
                               adjointVectorPosition = 0;}

  inline void Reset() {
    if (inputValues.size() != 0) {
      globalTape.reset();
      adjointVectorPosition = 0;
      inputValues.clear();
    }
  }

  inline void SetPreaccIn(const su2double &data) {
    if (PreaccActive) {
      if (data.isActive()) {
        localInputValues.push_back(data.getGradientData());
      }
    }
  }

  inline void SetPreaccIn(const su2double* data, const int size) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size; i++) {
        if (data[i].isActive()) {
          localInputValues.push_back(data[i].getGradientData());
        }
      }
    }
  }

  inline void SetPreaccIn(const su2double* const *data, const int size_x, const int size_y) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size_x; i++) {
        for (unsigned short j = 0; j < size_y; j++) {
          if (data[i][j].isActive()) {
            localInputValues.push_back(data[i][j].getGradientData());
          }
        }
      }
    }
  }

  inline void StartPreacc() {
    if (globalTape.isActive() && PreaccEnabled) {
      StartPosition = globalTape.getPosition();
      PreaccActive = true;
    }
  }

  inline void SetPreaccOut(su2double& data) {
    if (PreaccActive) {
      if (data.isActive()) {
        localOutputValues.push_back(&data);
      }
    }
  }

  inline void SetPreaccOut(su2double* data, const int size) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size; i++) {
        if (data[i].isActive()) {
          localOutputValues.push_back(&data[i]);
        }
      }
    }
  }

  inline void SetPreaccOut(su2double** data, const int size_x, const int size_y) {
    if (PreaccActive) {
      for (unsigned short i = 0; i < size_x; i++) {
        for (unsigned short j = 0; j < size_y; j++) {
          if (data[i][j].isActive()) {
            localOutputValues.push_back(&data[i][j]);
          }
        }
      }
    }
  }


  inline void delete_handler(void *handler) {
    CheckpointHandler *checkpoint = static_cast<CheckpointHandler*>(handler);
    checkpoint->clear();
  }
#else

  /*--- Default implementation if reverse mode is disabled ---*/

  inline void RegisterInput(su2double &data) {}

  inline void RegisterOutput(su2double& data) {}

  inline void StartRecording() {}

  inline void StopRecording() {}

  inline void ClearAdjoints() {}

  inline void ComputeAdjoint() {}

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
#endif
}
