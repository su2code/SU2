/*!
 * \file SU2_CFD.cpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
 * \author F. Palacios, T. Economon
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


#include "../include/SU2_CFD.hpp"

/* LIBXSMM include files, if supported. */
#ifdef HAVE_LIBXSMM
#include "libxsmm.h"
#endif

/* Include file, needed for the runtime NaN catching. */
//#include <fenv.h>

using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short nZone;
  char config_file_name[MAX_STRING_SIZE];
  bool fsi, turbo;
  bool dry_run = false;
  std::string filename = "default";
  
  /*--- Command line parsing ---*/
  
  CLI::App app{"SU2 v7.0.0 \"Blackbird\", The Open-Source CFD Code"};
  app.add_flag("-d,--dryrun", dry_run, "Enable dry run mode.\n" 
                                       "Only execute preprocessing steps using a dummy geometry.");
  app.add_option("configfile", filename, "A config file.")->check(CLI::ExistingFile);
  
  CLI11_PARSE(app, argc, argv)
  
  /*--- MPI initialization, and buffer setting ---*/
  
#ifdef HAVE_MPI
  int  buffsize;
  char *buffptr;
#ifdef HAVE_OMP
  int provided;
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
  
  CDriver *driver = NULL;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  strcpy(config_file_name, filename.c_str());

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation). ---*/

  CConfig *config = new CConfig(config_file_name, SU2_CFD);
  nZone  = config->GetnZone();
  fsi    = config->GetFSI_Simulation();
  turbo  = config->GetBoolTurbomachinery();

  /*--- First, given the basic information about the number of zones and the
   solver types from the config, instantiate the appropriate driver for the problem
   and perform all the preprocessing. ---*/
  
  if (!dry_run){
    
    if ((!config->GetMultizone_Problem() && (config->GetTime_Marching() != HARMONIC_BALANCE) && !turbo)
        || (turbo && config->GetDiscrete_Adjoint())) {
      
      /*--- Single zone problem: instantiate the single zone driver class. ---*/
      
      if (nZone > 1 ) {
        SU2_MPI::Error("The required solver doesn't support multizone simulations", CURRENT_FUNCTION);
      }
      if (config->GetDiscrete_Adjoint()) {

        driver = new CDiscAdjSinglezoneDriver(config_file_name, nZone, MPICommunicator);
      }
      else
        driver = new CSinglezoneDriver(config_file_name, nZone, MPICommunicator);
      
    }
    else if (config->GetMultizone_Problem() && !turbo && !fsi) {
      
    /*--- Multizone Drivers. ---*/

    if (config->GetDiscrete_Adjoint()) {

      driver = new CDiscAdjMultizoneDriver(config_file_name, nZone, MPICommunicator);

    }
    else {
      
      driver = new CMultizoneDriver(config_file_name, nZone, MPICommunicator);

    }
      
    } else if (config->GetTime_Marching() == HARMONIC_BALANCE) {
      
      /*--- Harmonic balance problem: instantiate the Harmonic Balance driver class. ---*/
      
      driver = new CHBDriver(config_file_name, nZone, MPICommunicator);
      
    } else if ((nZone == 2) && fsi) {
      
      bool stat_fsi = ((!config->GetTime_Domain()));
      bool disc_adj_fsi = (config->GetDiscrete_Adjoint());
      
      /*--- If the problem is a discrete adjoint FSI problem ---*/
      if (disc_adj_fsi) {
        if (stat_fsi) {
          driver = new CDiscAdjFSIDriver(config_file_name, nZone, MPICommunicator);
        }
        else {
          SU2_MPI::Error("WARNING: There is no discrete adjoint implementation for dynamic FSI. ", CURRENT_FUNCTION);
        }
      }
      
    } else {
      
      /*--- Multi-zone problem: instantiate the multi-zone driver class by default
            or a specialized driver class for a particular multi-physics problem. ---*/
      
      if (turbo) {
        
        driver = new CTurbomachineryDriver(config_file_name, nZone, MPICommunicator);
        
      } else {
        
        /*--- Instantiate the class for external aerodynamics ---*/
        
        driver = new CFluidDriver(config_file_name, nZone, MPICommunicator);
        
      }
      
    }
  } else {
    driver = new CDummyDriver(config_file_name, nZone, MPICommunicator);
  }
  
  delete config;
  config = NULL;
  
  /*--- Launch the main external loop of the solver ---*/

  driver->StartSolver();

  /*--- Postprocess all the containers, close history file, exit SU2 ---*/
  
  driver->Postprocessing();

  if (driver != NULL) delete driver;
  driver = NULL;
  
  /*---Finalize libxsmm, if supported. ---*/
#ifdef HAVE_LIBXSMM
  libxsmm_finalize();
#endif

  /*--- Finalize MPI parallelization ---*/
#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
