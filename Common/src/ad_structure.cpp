/*!
 * \file ad_structure.cpp
 * \brief Main subroutines for the algorithmic differentiation (AD) structure.
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

#include "../include/datatype_structure.hpp"

namespace AD {
#ifdef CODI_REVERSE_TYPE
  /*--- Initialization of the global variables ---*/

  int adjointVectorPosition = 0;

  std::vector<su2double::GradientData> inputValues;
  std::vector<su2double::GradientData> localInputValues;
  std::vector<su2double*> localOutputValues;

  su2double::TapeType& globalTape = su2double::getGlobalTape();
  su2double::TapeType::Position StartPosition, EndPosition;

  bool Status = false;
  bool PreaccActive = false;
  bool PreaccEnabled = true;

  void EndPreacc() {

    if(PreaccActive) {
      unsigned short iVarOut, iVarIn;
      unsigned short nVarOut, nVarIn;
      su2double::GradientData index_out, index_in;

      nVarOut = localOutputValues.size();
      nVarIn  = localInputValues.size();

      /*--- Store the current position of the tape ---*/

      EndPosition = globalTape.getPosition();

      /*--- Allocate local memory ---*/

      passivedouble* local_jacobi     = new passivedouble[nVarOut*nVarIn];
      unsigned short* nNonzero        = new unsigned short[nVarOut];

      /*--- Compute the local Jacobi matrix of the code between the start and end position
       * using the inputs and outputs declared with StartPreacc(...)/EndPreacc(...) ---*/

      for (iVarOut = 0; iVarOut < nVarOut; iVarOut++) {
        nNonzero[iVarOut] = 0;
        index_out = localOutputValues[iVarOut]->getGradientData();

        globalTape.setGradient(index_out, 1.0);
        globalTape.evaluate(EndPosition, StartPosition);

        for (iVarIn= 0; iVarIn < nVarIn; iVarIn++) {
          index_in =  localInputValues[iVarIn];
          local_jacobi[iVarOut*nVarIn+iVarIn] = globalTape.getGradient(index_in);
          if (local_jacobi[iVarOut*nVarIn+iVarIn] != 0.0) {
            nNonzero[iVarOut]++;
          }
          globalTape.setGradient(index_in, 0.0);
        }
        globalTape.setGradient(index_out, 0.0);
        globalTape.clearAdjoints(EndPosition, StartPosition);
      }

      /*--- Reset the tape to the starting position (to reuse the part of the tape) ---*/

      if (nVarOut > 0) {
        globalTape.reset(StartPosition);
      }

      /*--- For each output create a statement on the tape and push the corresponding Jacobi entries.
       * Note that the output variables need a new index since we did a reset of the tape section. ---*/

      for (iVarOut = 0; iVarOut < nVarOut; iVarOut++) {
        if (nNonzero[iVarOut] != 0){
          globalTape.store(localOutputValues[iVarOut]->getValue(), localOutputValues[iVarOut]->getGradientData(), nNonzero[iVarOut]);
          for (iVarIn = 0; iVarIn < nVarIn; iVarIn++) {
            index_in =  localInputValues[iVarIn];
           globalTape.pushJacobi(local_jacobi[iVarOut*nVarIn+iVarIn],
               local_jacobi[iVarOut*nVarIn+iVarIn], local_jacobi[iVarOut*nVarIn+iVarIn], index_in);
          }
        }
      }

      /*--- Clear local vectors and reset indicator ---*/


      localInputValues.clear();
      localOutputValues.clear();

      delete [] local_jacobi;
      delete [] nNonzero;

      PreaccActive = false;
    }
  }
#endif
}


/*--- If we compile under OSX we have to overload some of the operators for
 *   complex numbers to avoid the use of the standard operators
 *  (they use a lot of functions that are only defined for doubles) ---*/

#ifdef __APPLE__

namespace std{
  template<>
  su2double abs(const complex<su2double>& x){
    return sqrt(x.real()*x.real() + x.imag()*x.imag());
  }

  template<>
  complex<su2double> operator/(const complex<su2double>& x, const complex<su2double>& y){

    su2double d    = (y.real()*y.real() + y.imag()*y.imag());
    su2double real = (x.real()*y.real() + x.imag()*y.imag())/d;
    su2double imag = (x.imag()*y.real() - x.real()*y.imag())/d;

    return complex<su2double>(real, imag);

  }

  template<>
  complex<su2double> operator*(const complex<su2double>& x, const complex<su2double>& y){

    su2double real = (x.real()*y.real() - x.imag()*y.imag());
    su2double imag = (x.imag()*y.real() + x.real()*y.imag());

    return complex<su2double>(real, imag);

  }
}
#endif
