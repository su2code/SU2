/*!
 * \file CDummyDriver.cpp
 * \brief Dummy driver class for running the preprocessing without geometry preprocessing.
 * \author T. Albring
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/drivers/CDummyDriver.hpp"
#include "../../include/output/COutput.hpp"

CDummyDriver::CDummyDriver(char* confFile,
                         unsigned short val_nZone,
                         SU2_Comm MPICommunicator) : CDriver(confFile,
                                                             val_nZone,
                                                             MPICommunicator,
                                                             true) {
}

void CDummyDriver::StartSolver(){
  if (rank == MASTER_NODE){
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
    cout << endl;
    cout << "--------------------------------------------" << endl;
    cout << "No solver started. DRY_RUN option enabled. " << endl;
    cout << "--------------------------------------------" << endl;
  }

  for (iZone = 0; iZone < nZone; iZone++){
    output_container[iZone]->PrintVolumeFields();
    output_container[iZone]->PrintHistoryFields();
  }
}
