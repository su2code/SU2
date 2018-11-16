/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving multi-zone problems.
 * \author R. Sanchez, O. Burghardt
 * \version 6.0.1 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"

CDiscAdjMultizoneDriver::CDiscAdjMultizoneDriver(char* confFile,
                                      unsigned short val_nZone,
                                      unsigned short val_nDim,
                                      bool val_periodic,
                                      SU2_Comm MPICommunicator) : CMultizoneDriver(confFile,
                                                                                  val_nZone,
                                                                                  val_nDim,
                                                                                  val_periodic,
                                                                                  MPICommunicator) {

  RecordingState = NONE;
  unsigned short jZone;

  nInst = new unsigned short[nZone];

  direct_iteration = new CIteration**[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {

    nInst[iZone]            = 1;
    direct_iteration[iZone] = new CIteration*[nInst[iZone]];

    for(iInst = 0; iInst < nInst[iZone]; iInst++) {

      switch (config_container[iZone]->GetKind_Solver()) {

        case EULER: case NAVIER_STOKES: case RANS:
        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          direct_iteration[iZone][iInst] = new CFluidIteration(config_container[iZone]);
          break;

        default:
          SU2_MPI::Error("There is no discrete adjoint functionality for one of the specified solvers yet.", CURRENT_FUNCTION);
      }
    }
  }
}

CDiscAdjMultizoneDriver::~CDiscAdjMultizoneDriver(){

  for (iZone = 0; iZone < nZone; iZone++){
    for (iInst = 0; iInst < nInst[iInst]; iInst++){
      delete direct_iteration[iZone][iInst];
    }
    delete direct_iteration[iZone];
  }

  delete [] direct_iteration;
}
