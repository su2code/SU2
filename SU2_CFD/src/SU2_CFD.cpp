/*!
 * \file SU2_CFD.cpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../include/SU2_CFD.hpp"

/* LIBXSMM include files, if supported. */
#ifdef HAVE_LIBXSMM
#include "libxsmm.h"
#endif

/* Include file, needed for the runtime NaN catching. */
//#include <fenv.h>

using namespace std;

int main(int argc, char *argv[]) {

  char config_file_name[MAX_STRING_SIZE];
  bool dry_run = false;
  int num_threads = omp_get_max_threads();
  bool use_thread_mult = false;
  std::string filename = "default.cfg";

  /*--- Command line parsing ---*/

  CLI::App app{"SU2 v7.0.3 \"Blackbird\", The Open-Source CFD Code"};
  app.add_flag("-d,--dryrun", dry_run, "Enable dry run mode.\n"
                                       "Only execute preprocessing steps using a dummy geometry.");
  app.add_option("-t,--threads", num_threads, "Number of OpenMP threads per MPI rank.");
  app.add_flag("--thread_multiple", use_thread_mult, "Request MPI_THREAD_MULTIPLE thread support.");
  app.add_option("configfile", filename, "A config file.")->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv)

  omp_set_num_threads(num_threads);

  /*--- MPI initialization, and buffer setting ---*/

#ifdef HAVE_MPI
  int  buffsize;
  char *buffptr;
#ifdef HAVE_OMP
  int provided;
  if (use_thread_mult)
    SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  else
    SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  SU2_MPI::Init(&argc, &argv);
#endif
  SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
#else
  SU2_Comm MPICommunicator(0);
#endif

  /*--- Uncomment the following line if runtime NaN catching is desired. ---*/
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

  /*--- Initialize libxsmm, if supported. ---*/
#ifdef HAVE_LIBXSMM
  libxsmm_init();
#endif

  /*--- Create a pointer to the main SU2 Driver ---*/

  CDriver* driver = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  strcpy(config_file_name, filename.c_str());

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation). ---*/

  CConfig* config = new CConfig(config_file_name, SU2_CFD);
  unsigned short nZone = config->GetnZone();
  bool fsi = config->GetFSI_Simulation();
  bool turbo = config->GetBoolTurbomachinery();

  /*--- First, given the basic information about the number of zones and the
   solver types from the config, instantiate the appropriate driver for the problem
   and perform all the preprocessing. ---*/

  bool disc_adj = config->GetDiscrete_Adjoint();
  bool multizone = config->GetMultizone_Problem();
  bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);

  if (dry_run) {

    /*--- Dry Run. ---*/
    driver = new CDummyDriver(config_file_name, nZone, MPICommunicator);

  }
  else if ((!multizone && !harmonic_balance && !turbo) || (turbo && disc_adj)) {

    /*--- Generic single zone problem: instantiate the single zone driver class. ---*/
    if (nZone != 1 )
      SU2_MPI::Error("The required solver doesn't support multizone simulations", CURRENT_FUNCTION);

    if (disc_adj) {
      driver = new CDiscAdjSinglezoneDriver(config_file_name, nZone, MPICommunicator);
    }
    else {
      driver = new CSinglezoneDriver(config_file_name, nZone, MPICommunicator);
    }

  }
  else if (multizone && !turbo && !fsi) {

    /*--- Generic multizone problems. ---*/
    if (disc_adj) {
      driver = new CDiscAdjMultizoneDriver(config_file_name, nZone, MPICommunicator);
    }
    else {
      driver = new CMultizoneDriver(config_file_name, nZone, MPICommunicator);
    }

  }
  else if (harmonic_balance) {

    /*--- Harmonic balance problem: instantiate the Harmonic Balance driver class. ---*/
    driver = new CHBDriver(config_file_name, nZone, MPICommunicator);

  }
  else if (fsi && disc_adj) {

    /*--- Discrete adjoint FSI problem with the legacy driver. ---*/
    if (config->GetTime_Domain())
      SU2_MPI::Error("There is no discrete adjoint implementation for dynamic FSI. ", CURRENT_FUNCTION);

    if (nZone != 2)
      SU2_MPI::Error("The legacy discrete adjoint FSI driver only works for two-zone problems. ", CURRENT_FUNCTION);

    driver = new CDiscAdjFSIDriver(config_file_name, nZone, MPICommunicator);

  }
  else if (turbo) {

    /*--- Turbomachinery problem. ---*/
    driver = new CTurbomachineryDriver(config_file_name, nZone, MPICommunicator);

  }
  else {

    /*--- Instantiate the class for external aerodynamics by default. ---*/
    driver = new CFluidDriver(config_file_name, nZone, MPICommunicator);

  }

  delete config;
  config = nullptr;

  /*--- Launch the main external loop of the solver. ---*/

  driver->StartSolver();

  /*--- Postprocess all the containers, close history file, exit SU2. ---*/

  driver->Postprocessing();

  delete driver;
  driver = nullptr;

  /*---Finalize libxsmm, if supported. ---*/
#ifdef HAVE_LIBXSMM
  libxsmm_finalize();
#endif

  /*--- Finalize MPI parallelization. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif

  return EXIT_SUCCESS;

}
