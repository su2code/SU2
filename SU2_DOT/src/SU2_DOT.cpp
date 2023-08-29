/*!
 * \file SU2_DOT.cpp
 * \brief Main file of the Gradient Projection Code (SU2_DOT).
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../SU2_DEF/include/drivers/CDiscAdjDeformationDriver.hpp"

int main(int argc, char* argv[]) {
  char config_file_name[MAX_STRING_SIZE];

  /*--- MPI initialization. ---*/

#if defined(HAVE_OMP) && defined(HAVE_MPI)
  int provided;
  SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  SU2_MPI::Init(&argc, &argv);
#endif
  SU2_MPI::Comm comm = SU2_MPI::GetComm();

  /*--- Further initializations are placed in the constructor of CDriverBase, to ensure that they are also seen by the
   python wrapper. */

  /*--- Load in the number of zones and spatial dimensions in the mesh file
   (if no config file is specified, default.cfg is used). ---*/

  if (argc == 2) {
    strcpy(config_file_name, argv[1]);
  } else {
    strcpy(config_file_name, "default.cfg");
  }

  /*--- Initialize the mesh deformation driver. ---*/

  CDiscAdjDeformationDriver driver(config_file_name, comm);

  /*--- Preprocess the solver data. ---*/

  driver.Preprocess();

  /*--- Launch the main external loop of the solver. ---*/

  driver.Run();

  /*--- Postprocess all the containers, close history file, and exit SU2. ---*/

  driver.Finalize();

  /*--- Finalize AD. ---*/

  AD::Finalize();

  /*--- Finalize MPI parallelization. ---*/

  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}
