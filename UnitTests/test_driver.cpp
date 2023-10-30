/*!
 * \file test_driver.cpp
 * \brief The main entry point for unit tests (the main()).
 * \author C. Pederson
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

/*--- This macro tells Catch2 that we are defining a custom main (see
 * below) This should only be done in one file.  It should not be done
 * in any unit test files. ---*/
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "../../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../Common/include/option_structure.hpp"

int main(int argc, char* argv[]) {
  /*--- Startup MPI, if supported ---*/
#if defined(HAVE_OMP) && defined(HAVE_MPI)
  int provided;
  SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  SU2_MPI::Init(&argc, &argv);
#endif

  /*--- Run the test driver supplied by Catch ---*/
  int result = Catch::Session().run(argc, argv);

  /*--- Finalize MPI parallelization ---*/
  SU2_MPI::Finalize();

  return result;
}
