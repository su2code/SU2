/*!
 * \file CDriver.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez, F. Palacios
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/drivers/CDriver.hpp"
#include "../../include/definition_structure.hpp"

#include "../../../Common/include/geometry/CDummyGeometry.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../Common/include/geometry/CMultiGridGeometry.hpp"

#include "../../include/solvers/CSolverFactory.hpp"
#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"

#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"

#include "../../include/output/COutputLegacy.hpp"

#include "../../../Common/include/interface_interpolation/CInterpolator.hpp"
#include "../../../Common/include/interface_interpolation/CInterpolatorFactory.hpp"

#include "../../include/interfaces/cfd/CConservativeVarsInterface.hpp"
#include "../../include/interfaces/cfd/CMixingPlaneInterface.hpp"
#include "../../include/interfaces/cfd/CSlidingInterface.hpp"
#include "../../include/interfaces/cht/CConjugateHeatInterface.hpp"
#include "../../include/interfaces/fsi/CDisplacementsInterface.hpp"
#include "../../include/interfaces/fsi/CFlowTractionInterface.hpp"
#include "../../include/interfaces/fsi/CDiscAdjFlowTractionInterface.hpp"

#include "../../include/variables/CEulerVariable.hpp"
#include "../../include/variables/CIncEulerVariable.hpp"
#include "../../include/variables/CNEMOEulerVariable.hpp"

#include "../../include/numerics/template.hpp"
#include "../../include/numerics/transition.hpp"
#include "../../include/numerics/radiation.hpp"
#include "../../include/numerics/heat.hpp"
#include "../../include/numerics/flow/convection/roe.hpp"
#include "../../include/numerics/flow/convection/fds.hpp"
#include "../../include/numerics/flow/convection/fvs.hpp"
#include "../../include/numerics/flow/convection/cusp.hpp"
#include "../../include/numerics/flow/convection/hllc.hpp"
#include "../../include/numerics/flow/convection/ausm_slau.hpp"
#include "../../include/numerics/flow/convection/centered.hpp"
#include "../../include/numerics/flow/flow_diffusion.hpp"
#include "../../include/numerics/flow/flow_sources.hpp"
#include "../../include/numerics/NEMO/convection/roe.hpp"
#include "../../include/numerics/NEMO/convection/lax.hpp"
#include "../../include/numerics/NEMO/convection/ausm.hpp"
#include "../../include/numerics/NEMO/convection/ausmplusup2.hpp"
#include "../../include/numerics/NEMO/convection/ausmpwplus.hpp"
#include "../../include/numerics/NEMO/convection/msw.hpp"
#include "../../include/numerics/NEMO/NEMO_diffusion.hpp"
#include "../../include/numerics/NEMO/NEMO_sources.hpp"
#include "../../include/numerics/continuous_adjoint/adj_convection.hpp"
#include "../../include/numerics/continuous_adjoint/adj_diffusion.hpp"
#include "../../include/numerics/continuous_adjoint/adj_sources.hpp"
#include "../../include/numerics/scalar/scalar_convection.hpp"
#include "../../include/numerics/scalar/scalar_diffusion.hpp"
#include "../../include/numerics/scalar/scalar_sources.hpp"
#include "../../include/numerics/turbulent/turb_convection.hpp"
#include "../../include/numerics/turbulent/turb_diffusion.hpp"
#include "../../include/numerics/turbulent/turb_sources.hpp"
#include "../../include/numerics/species/species_convection.hpp"
#include "../../include/numerics/species/species_diffusion.hpp"
#include "../../include/numerics/species/species_sources.hpp"
#include "../../include/numerics/elasticity/CFEAElasticity.hpp"
#include "../../include/numerics/elasticity/CFEALinearElasticity.hpp"
#include "../../include/numerics/elasticity/CFEANonlinearElasticity.hpp"
#include "../../include/numerics/elasticity/nonlinear_models.hpp"

#include "../../include/integration/CIntegrationFactory.hpp"

#include "../../include/iteration/CIterationFactory.hpp"

#include "../../../Common/include/parallelization/omp_structure.hpp"

#include <cassert>

#ifdef VTUNEPROF
#include <ittnotify.h>
#endif
#include <fenv.h>

CDriver::CDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator, bool dummy_geo) :
  config_file_name(confFile), StartTime(0.0), StopTime(0.0), UsedTime(0.0),
  TimeIter(0), nZone(val_nZone), StopCalc(false), fsi(false), fem_solver(false), dry_run(dummy_geo) {

  /*--- Initialize Medipack (must also be here so it is initialized from python) ---*/
#ifdef HAVE_MPI
  #if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
    SU2_MPI::Init_AMPI();
  #endif
#endif

  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Start timer to track preprocessing for benchmarking. ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- Initialize containers with null --- */

  SetContainers_Null();

  /*--- Preprocessing of the config files. ---*/

  Input_Preprocessing(config_container, driver_config);

  /*--- Retrieve dimension from mesh file ---*/

  nDim = CConfig::GetnDim(config_container[ZONE_0]->GetMesh_FileName(),
                          config_container[ZONE_0]->GetMesh_FileFormat());

  /*--- Output preprocessing ---*/

  Output_Preprocessing(config_container, driver_config, output_container, driver_output);


  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone]    = new CGeometry**    [nInst[iZone]] ();
    iteration_container[iZone]   = new CIteration*    [nInst[iZone]] ();
    solver_container[iZone]      = new CSolver***     [nInst[iZone]] ();
    integration_container[iZone] = new CIntegration** [nInst[iZone]] ();
    numerics_container[iZone]    = new CNumerics****  [nInst[iZone]] ();
    grid_movement[iZone]         = new CVolumetricMovement* [nInst[iZone]] ();

    /*--- Allocate transfer and interpolation container --- */

    interface_container[iZone]    = new CInterface*[nZone] ();
    interpolator_container[iZone].resize(nZone);

    for (iInst = 0; iInst < nInst[iZone]; iInst++) {

      config_container[iZone]->SetiInst(iInst);

      /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
       based data structure is constructed, i.e. node and cell neighbors are
       identified and linked, face areas and volumes of the dual mesh cells are
       computed, and the multigrid levels are created using an agglomeration procedure. ---*/

      Geometrical_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], dry_run);

    }
  }

  /*--- Before we proceed with the zone loop we have to compute the wall distances.
     * This computation depends on all zones at once. ---*/
  if (rank == MASTER_NODE)
    cout << "Computing wall distances." << endl;

  CGeometry::ComputeWallDistance(config_container, geometry_container);

  for (iZone = 0; iZone < nZone; iZone++) {

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Definition of the solver class: solver_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS].
       The solver classes are specific to a particular set of governing equations,
       and they contain the subroutines with instructions for computing each spatial
       term of the PDE, i.e. loops over the edges to compute convective and viscous
       fluxes, loops over the nodes to compute source terms, and routines for
       imposing various boundary condition type for the PDE. ---*/

      Solver_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], solver_container[iZone][iInst]);

      /*--- Definition of the numerical method class:
       numerics_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
       The numerics class contains the implementation of the numerical methods for
       evaluating convective or viscous fluxes between any two nodes in the edge-based
       data structure (centered, upwind, galerkin), as well as any source terms
       (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/

      Numerics_Preprocessing(config_container[iZone], geometry_container[iZone][iInst],
                             solver_container[iZone][iInst], numerics_container[iZone][iInst]);

      /*--- Definition of the integration class: integration_container[#ZONES][#INSTANCES][#EQ_SYSTEMS].
       The integration class orchestrates the execution of the spatial integration
       subroutines contained in the solver class (including multigrid) for computing
       the residual at each node, R(U) and then integrates the equations to a
       steady state or time-accurately. ---*/

      Integration_Preprocessing(config_container[iZone], solver_container[iZone][iInst][MESH_0],
                                integration_container[iZone][iInst]);

      /*--- Instantiate the type of physics iteration to be executed within each zone. For
       example, one can execute the same physics across multiple zones (mixing plane),
       different physics in different zones (fluid-structure interaction), or couple multiple
       systems tightly within a single zone by creating a new iteration class (e.g., RANS). ---*/

      Iteration_Preprocessing(config_container[iZone], iteration_container[iZone][iInst]);

      /*--- Dynamic mesh processing.  ---*/

      DynamicMesh_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], solver_container[iZone][iInst],
                                iteration_container[iZone][iInst], grid_movement[iZone][iInst], surface_movement[iZone]);

      /*--- Static mesh processing.  ---*/

      StaticMesh_Preprocessing(config_container[iZone], geometry_container[iZone][iInst]);

    }

  }

  /*! --- Compute the wall distance again to correctly compute the derivatives if we are running direct diff mode --- */
  if (driver_config->GetDirectDiff() == D_DESIGN){
    CGeometry::ComputeWallDistance(config_container, geometry_container);
  }

  /*--- Definition of the interface and transfer conditions between different zones. ---*/

  if (nZone > 1) {
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Multizone Interface Preprocessing -------------------" << endl;

    Interface_Preprocessing(config_container, solver_container, geometry_container,
                            interface_types, interface_container, interpolator_container);
  }

  if (fsi) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        Solver_Restart(solver_container[iZone][iInst], geometry_container[iZone][iInst],
                       config_container[iZone], true);
      }
    }
  }

  if (config_container[ZONE_0]->GetBoolTurbomachinery()){
    if (rank == MASTER_NODE)
      cout << endl <<"---------------------- Turbomachinery Preprocessing ---------------------" << endl;

    Turbomachinery_Preprocessing(config_container, geometry_container, solver_container, interface_container);
  }


  PythonInterface_Preprocessing(config_container, geometry_container, solver_container);


  /*--- Preprocessing time is reported now, but not included in the next compute portion. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc    = UsedTime;
  UsedTimeCompute    = 0.0;
  UsedTimeOutput     = 0.0;
  IterCount          = 0;
  OutputCount        = 0;
  MDOFs              = 0.0;
  MDOFsDomain        = 0.0;
  Mpoints            = 0.0;
  MpointsDomain      = 0.0;
  for (iZone = 0; iZone < nZone; iZone++) {
    Mpoints       += geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPoint()/(1.0e6);
    MpointsDomain += geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPointDomain()/(1.0e6);
    MDOFs         += DOFsPerPoint*geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPoint()/(1.0e6);
    MDOFsDomain   += DOFsPerPoint*geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPointDomain()/(1.0e6);
  }

  /*--- Reset timer for compute/output performance benchmarking. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc = UsedTime;

  /*--- Reset timer for compute performance benchmarking. ---*/

  StartTime = SU2_MPI::Wtime();

}

void CDriver::SetContainers_Null(){

  /*--- Create pointers to all of the classes that may be used throughout
   the SU2_CFD code. In general, the pointers are instantiated down a
   hierarchy over all zones, multigrid levels, equation sets, and equation
   terms as described in the comments below. ---*/

  ConvHist_file                  = nullptr;
  iteration_container            = nullptr;
  output_container               = nullptr;
  integration_container          = nullptr;
  geometry_container             = nullptr;
  solver_container               = nullptr;
  numerics_container             = nullptr;
  config_container               = nullptr;
  surface_movement               = nullptr;
  grid_movement                  = nullptr;
  FFDBox                         = nullptr;
  interface_container            = nullptr;
  interface_types                = nullptr;
  nInst                          = nullptr;

  /*--- Definition and of the containers for all possible zones. ---*/

  iteration_container            = new CIteration**[nZone] ();
  solver_container               = new CSolver****[nZone] ();
  integration_container          = new CIntegration***[nZone] ();
  numerics_container             = new CNumerics*****[nZone] ();
  config_container               = new CConfig*[nZone] ();
  geometry_container             = new CGeometry***[nZone] ();
  surface_movement               = new CSurfaceMovement*[nZone] ();
  grid_movement                  = new CVolumetricMovement**[nZone] ();
  FFDBox                         = new CFreeFormDefBox**[nZone] ();
  interpolator_container.resize(nZone);
  interface_container            = new CInterface**[nZone] ();
  interface_types                = new unsigned short*[nZone] ();
  output_container               = new COutput*[nZone] ();
  nInst                          = new unsigned short[nZone] ();
  driver_config                  = nullptr;
  driver_output                  = nullptr;

  for (iZone = 0; iZone < nZone; iZone++) {
    interface_types[iZone] = new unsigned short[nZone];
    nInst[iZone] = 1;
  }

  strcpy(runtime_file_name, "runtime.dat");

}


void CDriver::Postprocessing() {

  const bool wrt_perf = config_container[ZONE_0]->GetWrt_Performance();

    /*--- Output some information to the console. ---*/

  if (rank == MASTER_NODE) {

    /*--- Print out the number of non-physical points and reconstructions ---*/

    if (config_container[ZONE_0]->GetNonphysical_Points() > 0)
      cout << "Warning: there are " << config_container[ZONE_0]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
    if (config_container[ZONE_0]->GetNonphysical_Reconstr() > 0)
      cout << "Warning: " << config_container[ZONE_0]->GetNonphysical_Reconstr() << " reconstructed states for upwinding are non-physical." << endl;
  }

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Numerics_Postprocessing(numerics_container[iZone], solver_container[iZone][iInst],
          geometry_container[iZone][iInst], config_container[iZone], iInst);
    }
    delete [] numerics_container[iZone];
  }
  delete [] numerics_container;
  if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Integration_Postprocessing(integration_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    }
    delete [] integration_container[iZone];
  }
  delete [] integration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIntegration container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Solver_Postprocessing(solver_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    }
    delete [] solver_container[iZone];
  }
  delete [] solver_container;
  if (rank == MASTER_NODE) cout << "Deleted CSolver container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      delete iteration_container[iZone][iInst];
    delete [] iteration_container[iZone];
  }
  delete [] iteration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIteration container." << endl;

  if (interface_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (interface_container[iZone] != nullptr) {
        for (unsigned short jZone = 0; jZone < nZone; jZone++)
          delete interface_container[iZone][jZone];
        delete [] interface_container[iZone];
      }
    }
    delete [] interface_container;
    if (rank == MASTER_NODE) cout << "Deleted CInterface container." << endl;
  }

  if (interface_types != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++)
      delete [] interface_types[iZone];
    delete [] interface_types;
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone] != nullptr) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        for (unsigned short iMGlevel = 0; iMGlevel < config_container[iZone]->GetnMGLevels()+1; iMGlevel++)
          delete geometry_container[iZone][iInst][iMGlevel];
        delete [] geometry_container[iZone][iInst];
      }
      delete [] geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete [] FFDBox[iZone];
  }
  delete [] FFDBox;
  if (rank == MASTER_NODE) cout << "Deleted CFreeFormDefBox class." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete surface_movement[iZone];
  }
  delete [] surface_movement;
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      delete grid_movement[iZone][iInst];
    delete [] grid_movement[iZone];
  }
  delete [] grid_movement;
  if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;

  /*--- Output profiling information ---*/
  // Note that for now this is called only by a single thread, but all
  // necessary variables have been made thread private for safety (tick/tock)!!

  config_container[ZONE_0]->SetProfilingCSV();
  config_container[ZONE_0]->GEMMProfilingCSV();

  /*--- Deallocate config container ---*/
  if (config_container!= nullptr) {
    for (iZone = 0; iZone < nZone; iZone++)
      delete config_container[iZone];
    delete [] config_container;
  }
  delete driver_config;
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  delete [] nInst;
  if (rank == MASTER_NODE) cout << "Deleted nInst container." << endl;

  /*--- Deallocate output container ---*/

  if (output_container!= nullptr) {
    for (iZone = 0; iZone < nZone; iZone++)
      delete output_container[iZone];
    delete [] output_container;
  }

  delete driver_output;

  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl;


  /*--- Stop the timer and output the final performance summary. ---*/

  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime-StartTime;
  UsedTimeCompute += UsedTime;

  if ((rank == MASTER_NODE) && (wrt_perf)) {
    su2double TotalTime = UsedTimePreproc + UsedTimeCompute + UsedTimeOutput;
    cout.precision(6);
    cout << endl << endl <<"-------------------------- Performance Summary --------------------------" << endl;
    cout << "Simulation totals:" << endl;
    cout << setw(25) << "Wall-clock time (hrs):" << setw(12) << (TotalTime)/(60.0*60.0) << " | ";
    cout << setw(20) << "Core-hrs:" << setw(12) << size*TotalTime/(60.0*60.0) << endl;
    cout << setw(25) << "Cores:" << setw(12) << size << " | ";
    cout << setw(20) << "DOFs/point:" << setw(12) << DOFsPerPoint << endl;
    cout << setw(25) << "Points/core:" << setw(12) << 1.0e6*MpointsDomain/size << " | ";
    cout << setw(20) << "Ghost points/core:" << setw(12) << 1.0e6*(Mpoints-MpointsDomain)/size << endl;
    cout << setw(25) << "Ghost/Owned Point Ratio:" << setw(12) << (Mpoints-MpointsDomain)/MpointsDomain << " | " << endl;
    cout << endl;
    cout << "Preprocessing phase:" << endl;
    cout << setw(25) << "Preproc. Time (s):"  << setw(12)<< UsedTimePreproc << " | ";
    cout << setw(20) << "Preproc. Time (%):" << setw(12)<< ((UsedTimePreproc * 100.0) / (TotalTime)) << endl;
    cout << endl;
    cout << "Compute phase:" << endl;
    cout << setw(25) << "Compute Time (s):"  << setw(12)<< UsedTimeCompute << " | ";
    cout << setw(20) << "Compute Time (%):" << setw(12)<< ((UsedTimeCompute * 100.0) / (TotalTime)) << endl;
    cout << setw(25) << "Iteration count:"  << setw(12)<< IterCount << " | ";
    if (IterCount != 0) {
      cout << setw(20) << "Avg. s/iter:" << setw(12)<< UsedTimeCompute/IterCount << endl;
      cout << setw(25) << "Core-s/iter/Mpoints:" << setw(12)<< size*UsedTimeCompute/IterCount/Mpoints << " | ";
      cout << setw(20) << "Mpoints/s:" << setw(12)<< Mpoints*IterCount/UsedTimeCompute << endl;
    } else cout << endl;
    cout << endl;
    cout << "Output phase:" << endl;
    cout << setw(25) << "Output Time (s):"  << setw(12)<< UsedTimeOutput << " | ";
    cout << setw(20) << "Output Time (%):" << setw(12)<< ((UsedTimeOutput * 100.0) / (TotalTime)) << endl;
    cout << setw(25) << "Output count:" << setw(12)<< OutputCount << " | ";
    if (OutputCount != 0) {
      cout << setw(20)<< "Avg. s/output:" << setw(12)<< UsedTimeOutput/OutputCount << endl;
      if (BandwidthSum > 0) {
        cout << setw(25)<< "Restart Aggr. BW (MB/s):" << setw(12)<< BandwidthSum/OutputCount << " | ";
        cout << setw(20)<< "MB/s/core:" << setw(12)<< BandwidthSum/OutputCount/size << endl;
      }
    } else cout << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_CFD) ------------------------" << endl << endl;

}


void CDriver::Input_Preprocessing(CConfig **&config, CConfig *&driver_config) {

  char zone_file_name[MAX_STRING_SIZE];

  /*--- Initialize the configuration of the driver ---*/

  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_CFD, false);

  for (iZone = 0; iZone < nZone; iZone++) {

    if (rank == MASTER_NODE){
      cout  << endl << "Parsing config file for zone " << iZone << endl;
    }
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (driver_config->GetnConfigFiles() > 0){

      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_CFD, iZone, nZone, true);
    }
    else{
      config[iZone] = new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_CFD, iZone, nZone, true);
    }

    /*--- Set the MPI communicator ---*/

    config[iZone]->SetMPICommunicator(SU2_MPI::GetComm());
  }


  /*--- Set the multizone part of the problem. ---*/
  if (driver_config->GetMultizone_Problem()){
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }

  /*--- Determine whether or not the FEM solver is used, which decides the type of
   *    geometry classes that are instantiated. Only adapted for single-zone problems ---*/

  fem_solver = config_container[ZONE_0]->GetFEMSolver();

  fsi = config_container[ZONE_0]->GetFSI_Simulation();
}

void CDriver::Geometrical_Preprocessing(CConfig* config, CGeometry **&geometry, bool dummy){

  if (!dummy){
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Geometry Preprocessing ( Zone " << config->GetiZone() <<" ) -------------------" << endl;

    if( fem_solver ) {
      switch( config->GetKind_FEM_Flow() ) {
        case DG: {
            Geometrical_Preprocessing_DGFEM(config, geometry);
            break;
          }
      }
    }
    else {
      Geometrical_Preprocessing_FVM(config, geometry);
    }
  } else {
    if (rank == MASTER_NODE)
      cout << endl <<"-------------------------- Using Dummy Geometry -------------------------" << endl;

    unsigned short iMGlevel;

    geometry = new CGeometry*[config->GetnMGLevels()+1] ();

    if (!fem_solver){
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        geometry[iMGlevel] = new CDummyGeometry(config);
      }
    } else {
      geometry[MESH_0] = new CDummyMeshFEM_DG(config);
    }

    nDim = geometry[MESH_0]->GetnDim();
  }

  /*--- Computation of positive surface area in the z-plane which is used for
     the calculation of force coefficient (non-dimensionalization). ---*/

  geometry[MESH_0]->SetPositive_ZArea(config);

  /*--- Set the actuator disk boundary conditions, if necessary. ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    geometry[iMesh]->MatchActuator_Disk(config);
  }

  /*--- If we have any periodic markers in this calculation, we must
       match the periodic points found on both sides of the periodic BC.
       Note that the current implementation requires a 1-to-1 matching of
       periodic points on the pair of periodic faces after the translation
       or rotation is taken into account. ---*/

  if ((config->GetnMarker_Periodic() != 0) && !fem_solver) {
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

      /*--- Note that we loop over pairs of periodic markers individually
           so that repeated nodes on adjacent periodic faces are properly
           accounted for in multiple places. ---*/

      for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
        geometry[iMesh]->MatchPeriodic(config, iPeriodic);
      }

      /*--- For Streamwise Periodic flow, find a unique reference node on the dedicated inlet marker. ---*/
      if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)
        geometry[iMesh]->FindUniqueNode_PeriodicBound(config);

      /*--- Initialize the communication framework for the periodic BCs. ---*/
      geometry[iMesh]->PreprocessPeriodicComms(geometry[iMesh], config);

    }
  }

  /*--- If activated by the compile directive, perform a partition analysis. ---*/
#if PARTITION
  if (!dummy){
    if( fem_solver ) Partition_Analysis_FEM(geometry[MESH_0], config);
    else Partition_Analysis(geometry[MESH_0], config);
  }
#endif

  /*--- Check if Euler & Symmetry markers are straight/plane. This information
        is used in the Euler & Symmetry boundary routines. ---*/
  if((config_container[iZone]->GetnMarker_Euler() != 0 ||
     config_container[iZone]->GetnMarker_SymWall() != 0) &&
     !fem_solver) {

    if (rank == MASTER_NODE)
      cout << "Checking if Euler & Symmetry markers are straight/plane:" << endl;

    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++)
      geometry_container[iZone][iInst][iMesh]->ComputeSurf_Straightness(config_container[iZone], (iMesh==MESH_0) );

  }

}

void CDriver::Geometrical_Preprocessing_FVM(CConfig *config, CGeometry **&geometry) {

  unsigned short iZone = config->GetiZone(), iMGlevel;
  unsigned short requestedMGlevels = config->GetnMGLevels();
  const bool fea = config->GetStructuralProblem();

  /*--- Definition of the geometry class to store the primal grid in the partitioning process.
   *    All ranks process the grid and call ParMETIS for partitioning ---*/

  CGeometry *geometry_aux = new CPhysicalGeometry(config, iZone, nZone);

  /*--- Set the dimension --- */

  nDim = geometry_aux->GetnDim();

  /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

  geometry_aux->SetColorGrid_Parallel(config);

  /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

  geometry = new CGeometry *[config->GetnMGLevels()+1] ();

  /*--- Build the grid data structures using the ParMETIS coloring. ---*/

  geometry[MESH_0] = new CPhysicalGeometry(geometry_aux, config);

  /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

  delete geometry_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetBoundaries(config);

  /*--- Compute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();

  /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/

  if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  geometry[MESH_0]->SetRCM_Ordering(config);

  /*--- recompute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();

  /*--- Compute elements surrounding elements ---*/

  if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetElement_Connectivity();

  /*--- Check the orientation before computing geometrical quantities ---*/

  geometry[MESH_0]->SetBoundVolume();
  if (config->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[MESH_0]->Check_IntElem_Orientation(config);
    geometry[MESH_0]->Check_BoundElem_Orientation(config);
  }

  /*--- Create the edge structure ---*/

  if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);

  /*--- Create the control volume structures ---*/

  if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
  SU2_OMP_PARALLEL {
    geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
    geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);
  }
  END_SU2_OMP_PARALLEL

  /*--- Visualize a dual control volume if requested ---*/

  if ((config->GetVisualize_CV() >= 0) &&
      (config->GetVisualize_CV() < (long)geometry[MESH_0]->GetGlobal_nPointDomain()))
    geometry[MESH_0]->VisualizeControlVolume(config);

  /*--- Identify closest normal neighbor ---*/

  if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
  geometry[MESH_0]->FindNormal_Neighbor(config);

  /*--- Store the global to local mapping. ---*/

  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry[MESH_0]->SetGlobal_to_Local_Point();

  /*--- Compute the surface curvature ---*/

  if (!fea) {
    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry[MESH_0]->ComputeSurf_Curvature(config);
  }

  /*--- Compute the global surface areas for all markers. ---*/

  geometry[MESH_0]->ComputeSurfaceAreaCfgFile(config);

  /*--- Check for periodicity and disable MG if necessary. ---*/

  if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
  geometry[MESH_0]->Check_Periodicity(config);

  /*--- Compute mesh quality statistics on the fine grid. ---*/

  if (!fea) {
    if (rank == MASTER_NODE)
      cout << "Computing mesh quality statistics for the dual control volumes." << endl;
    geometry[MESH_0]->ComputeMeshQualityStatistics(config);
  }

  geometry[MESH_0]->SetMGLevel(MESH_0);
  if ((config->GetnMGLevels() != 0) && (rank == MASTER_NODE))
    cout << "Setting the multigrid structure." << endl;

  /*--- Loop over all the new grid ---*/

  for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Create main agglomeration structure ---*/

    geometry[iMGlevel] = new CMultiGridGeometry(geometry[iMGlevel-1], config, iMGlevel);

    /*--- Compute points surrounding points. ---*/

    geometry[iMGlevel]->SetPoint_Connectivity(geometry[iMGlevel-1]);

    /*--- Create the edge structure ---*/

    geometry[iMGlevel]->SetEdges();
    geometry[iMGlevel]->SetVertex(geometry[iMGlevel-1], config);

    /*--- Create the control volume structures ---*/

    geometry[iMGlevel]->SetControlVolume(geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetBoundControlVolume(geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGlevel-1]);

    /*--- Find closest neighbor to a surface point ---*/

    geometry[iMGlevel]->FindNormal_Neighbor(config);

    /*--- Store our multigrid index. ---*/

    geometry[iMGlevel]->SetMGLevel(iMGlevel);

    /*--- Protect against the situation that we were not able to complete
       the agglomeration for this level, i.e., there weren't enough points.
       We need to check if we changed the total number of levels and delete
       the incomplete CMultiGridGeometry object. ---*/

    if (config->GetnMGLevels() != requestedMGlevels) {
      delete geometry[iMGlevel];
      geometry[iMGlevel] = nullptr;
      break;
    }

  }

  if (config->GetWrt_MultiGrid()) geometry[MESH_0]->ColorMGLevels(config->GetnMGLevels(), geometry);

  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/

  if ((config->GetTime_Marching() != TIME_MARCHING::STEADY) && config->GetDynamic_Grid()) {
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      /*--- Update cell volume ---*/
      geometry[iMGlevel]->nodes->SetVolume_n();
      geometry[iMGlevel]->nodes->SetVolume_nM1();

      if (config->GetGrid_Movement()) {
        /*--- Update point coordinates ---*/
        geometry[iMGlevel]->nodes->SetCoord_n();
        geometry[iMGlevel]->nodes->SetCoord_n1();
      }
    }
  }


  /*--- Create the data structure for MPI point-to-point communications. ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
    geometry[iMGlevel]->PreprocessP2PComms(geometry[iMGlevel], config);


  /*--- Perform a few preprocessing routines and communications. ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Compute the max length. ---*/

    if (!fea) {
      if ((rank == MASTER_NODE) && (iMGlevel == MESH_0))
        cout << "Finding max control volume width." << endl;
      geometry[iMGlevel]->SetMaxLength(config);
    }

    /*--- Communicate the number of neighbors. This is needed for
         some centered schemes and for multigrid in parallel. ---*/

    if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (iMGlevel == MESH_0))
      cout << "Communicating number of neighbors." << endl;
    geometry[iMGlevel]->InitiateComms(geometry[iMGlevel], config, NEIGHBORS);
    geometry[iMGlevel]->CompleteComms(geometry[iMGlevel], config, NEIGHBORS);
  }

}

void CDriver::Geometrical_Preprocessing_DGFEM(CConfig* config, CGeometry **&geometry) {

  /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
  /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

  CGeometry *geometry_aux = new CPhysicalGeometry(config, iZone, nZone);

  /*--- Set the dimension --- */

  nDim = geometry_aux->GetnDim();

  /*--- For the FEM solver with time-accurate local time-stepping, use
       a dummy solver class to retrieve the initial flow state. ---*/

  CSolver *solver_aux = new CFEM_DG_EulerSolver(config, nDim, MESH_0);

  /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

  geometry_aux->SetColorFEMGrid_Parallel(config);

  /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

  geometry = new CGeometry *[config->GetnMGLevels()+1] ();

  geometry[MESH_0] = new CMeshFEM_DG(geometry_aux, config);

  /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

  delete geometry_aux;
  delete solver_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetBoundaries(config);

  /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
       define all virtual functions in the base class CGeometry. ---*/
  CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry[MESH_0]);

  /*--- Determine the standard elements for the volume elements. ---*/
  if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
  DGMesh->CreateStandardVolumeElements(config);

  /*--- Create the face information needed to compute the contour integral
       for the elements in the Discontinuous Galerkin formulation. ---*/
  if (rank == MASTER_NODE) cout << "Creating face information." << endl;
  DGMesh->CreateFaces(config);

  /*--- Compute the metric terms of the volume elements. ---*/
  if (rank == MASTER_NODE) cout << "Computing metric terms volume elements." << endl;
  DGMesh->MetricTermsVolumeElements(config);

  /*--- Compute the metric terms of the surface elements. ---*/
  if (rank == MASTER_NODE) cout << "Computing metric terms surface elements." << endl;
  DGMesh->MetricTermsSurfaceElements(config);

  /*--- Compute a length scale of the volume elements. ---*/
  if (rank == MASTER_NODE) cout << "Computing length scale volume elements." << endl;
  DGMesh->LengthScaleVolumeElements();

  /*--- Compute the coordinates of the integration points. ---*/
  if (rank == MASTER_NODE) cout << "Computing coordinates of the integration points." << endl;
  DGMesh->CoordinatesIntegrationPoints();

  /*--- Compute the coordinates of the location of the solution DOFs. This is different
            from the grid points when a different polynomial degree is used to represent the
            geometry and solution. ---*/
  if (rank == MASTER_NODE) cout << "Computing coordinates of the solution DOFs." << endl;
  DGMesh->CoordinatesSolDOFs();

  /*--- Perform the preprocessing tasks when wall functions are used. ---*/
  if (rank == MASTER_NODE) cout << "Preprocessing for the wall functions. " << endl;
  DGMesh->WallFunctionPreprocessing(config);

  /*--- Store the global to local mapping. ---*/
  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local DOF index." << endl;
  geometry[MESH_0]->SetGlobal_to_Local_Point();


  /*--- Loop to create the coarser grid levels. ---*/

  for(unsigned short iMGlevel=1; iMGlevel<=config->GetnMGLevels(); iMGlevel++) {

    SU2_MPI::Error("Geometrical_Preprocessing_DGFEM: Coarse grid levels not implemented yet.",
                   CURRENT_FUNCTION);
  }

}

void CDriver::Solver_Preprocessing(CConfig* config, CGeometry** geometry, CSolver ***&solver) {

  MAIN_SOLVER kindSolver = config->GetKind_Solver();

  if (rank == MASTER_NODE)
    cout << endl <<"-------------------- Solver Preprocessing ( Zone " << config->GetiZone() <<" ) --------------------" << endl;

  solver = new CSolver**[config->GetnMGLevels()+1] ();

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++){
    solver[iMesh] = CSolverFactory::CreateSolverContainer(kindSolver, config, geometry[iMesh], iMesh);
  }

  /*--- Count the number of DOFs per solution point. ---*/

  DOFsPerPoint = 0;

  for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++)
    if (solver[MESH_0][iSol]) DOFsPerPoint += solver[MESH_0][iSol]->GetnVar();

  /*--- Restart solvers, for FSI the geometry cannot be updated because the interpolation classes
   * should always use the undeformed mesh (otherwise the results would not be repeatable). ---*/

  if (!fsi) Solver_Restart(solver, geometry, config, true);

  /*--- Set up any necessary inlet profiles ---*/

  Inlet_Preprocessing(solver, geometry, config);

}

void CDriver::Inlet_Preprocessing(CSolver ***solver, CGeometry **geometry,
                                  CConfig *config) const {

  /*--- Adjust iteration number for unsteady restarts. ---*/

  const bool adjoint = config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint();

  int val_iter = 0;

  if (config->GetTime_Domain()) {
    val_iter = adjoint? config->GetUnst_AdjointIter() : config->GetRestart_Iter();
    val_iter -= 1;
    if (!adjoint && config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
      val_iter -= 1;
    if (!adjoint && !config->GetRestart()) val_iter = 0;
  }

  /*--- Load inlet profile files for any of the active solver containers.
   Note that these routines fill the fine grid data structures for the markers
   and restrict values down to all coarser MG levels. ---*/

  if (config->GetInlet_Profile_From_File()) {

    /*--- Use LoadInletProfile() routines for the particular solver. ---*/

    if (rank == MASTER_NODE) {
      cout << endl;
      cout << "Reading inlet profile from file: ";
      cout << config->GetInlet_FileName() << endl;
    }

    if (solver[MESH_0][FLOW_SOL]) {
      solver[MESH_0][FLOW_SOL]->LoadInletProfile(geometry, solver, config, val_iter, FLOW_SOL, INLET_FLOW);
    }
    if (solver[MESH_0][TURB_SOL]) {
      solver[MESH_0][TURB_SOL]->LoadInletProfile(geometry, solver, config, val_iter, TURB_SOL, INLET_FLOW);
    }
    if (solver[MESH_0][SPECIES_SOL]) {
      solver[MESH_0][SPECIES_SOL]->LoadInletProfile(geometry, solver, config, val_iter, SPECIES_SOL, INLET_FLOW);
    }

    /*--- Exit if profiles were requested for a solver that is not available. ---*/

    if (!config->GetFluidProblem()) {
      SU2_MPI::Error(string("Inlet profile specification via file (C++) has not been \n") +
                     string("implemented yet for this solver.\n") +
                     string("Please set SPECIFIED_INLET_PROFILE= NO and try again."), CURRENT_FUNCTION);
    }

  } else {

    /*--- Uniform inlets or python-customized inlets ---*/

    /* --- Initialize quantities for inlet boundary
     * This routine does not check if the python wrapper is being used to
     * set custom boundary conditions.  This is intentional; the
     * default values for python custom BCs are initialized with the default
     * values specified in the config (avoiding non physical values) --- */

    for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (solver[iMesh][FLOW_SOL]) solver[iMesh][FLOW_SOL]->SetUniformInlet(config, iMarker);
        if (solver[iMesh][TURB_SOL]) solver[iMesh][TURB_SOL]->SetUniformInlet(config, iMarker);
        if (solver[iMesh][SPECIES_SOL]) solver[iMesh][SPECIES_SOL]->SetUniformInlet(config, iMarker);
      }
    }

  }

}

void CDriver::Solver_Restart(CSolver ***solver, CGeometry **geometry,
                             CConfig *config, bool update_geo) {

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  const bool restart = config->GetRestart();
  const bool restart_flow = config->GetRestart_Flow();

  /*--- Adjust iteration number for unsteady restarts. ---*/

  int val_iter = 0;

  const bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());
  const bool time_domain = config->GetTime_Domain();
  const bool dt_step_2nd = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) &&
                           !config->GetStructuralProblem() && !config->GetFEMSolver() &&
                           !adjoint && time_domain;

  if (time_domain) {
    if (adjoint) val_iter = config->GetUnst_AdjointIter() - 1;
    else val_iter = config->GetRestart_Iter() - 1 - dt_step_2nd;
  }

  /*--- Restart direct solvers. ---*/

  if (restart || restart_flow) {
    for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
      auto sol = solver[MESH_0][iSol];
      if (sol && !sol->GetAdjoint()) {
        /*--- Note that the mesh solver always loads the most recent file (and not -2). ---*/
        SU2_OMP_PARALLEL_(if(sol->GetHasHybridParallel()))
        sol->LoadRestart(geometry, solver, config, val_iter + (iSol==MESH_SOL && dt_step_2nd), update_geo);
        END_SU2_OMP_PARALLEL
      }
    }
  }

  /*--- Restart adjoint solvers. ---*/

  if (restart) {
    if ((config->GetKind_Solver() == MAIN_SOLVER::TEMPLATE_SOLVER) ||
        (config->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS && !config->GetFrozen_Visc_Cont())) {
      SU2_MPI::Error("A restart capability has not been implemented yet for this solver.\n"
                     "Please set RESTART_SOL= NO and try again.", CURRENT_FUNCTION);
    }

    for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
      auto sol = solver[MESH_0][iSol];
      if (sol && sol->GetAdjoint())
        sol->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
  }

}

void CDriver::Solver_Postprocessing(CSolver ****solver, CGeometry **geometry,
                                    CConfig *config, unsigned short val_iInst) {

  for (int iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++){
      delete solver[val_iInst][iMGlevel][iSol];
    }
    delete [] solver[val_iInst][iMGlevel];
  }
  delete [] solver[val_iInst];

  CSolverFactory::ClearSolverMeta();

}

void CDriver::Integration_Preprocessing(CConfig *config, CSolver **solver, CIntegration **&integration) const {

  if (rank == MASTER_NODE)
    cout << endl <<"----------------- Integration Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  MAIN_SOLVER kindMainSolver = config->GetKind_Solver();

  integration = CIntegrationFactory::CreateIntegrationContainer(kindMainSolver, solver);

}

void CDriver::Integration_Postprocessing(CIntegration ***integration, CGeometry **geometry, CConfig *config, unsigned short val_iInst) {

  for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++){
    delete integration[val_iInst][iSol];
  }

  delete [] integration[val_iInst];

}

template <class Indices>
void CDriver::InstantiateTurbulentNumerics(unsigned short nVar_Turb, int offset, const CConfig *config,
                                           const CSolver* turb_solver, CNumerics ****&numerics) const {
  const int conv_term = CONV_TERM + offset;
  const int visc_term = VISC_TERM + offset;

  const int source_first_term = SOURCE_FIRST_TERM + offset;
  const int source_second_term = SOURCE_SECOND_TERM + offset;

  const int conv_bound_term = CONV_BOUND_TERM + offset;
  const int visc_bound_term = VISC_BOUND_TERM + offset;

  bool spalart_allmaras, neg_spalart_allmaras, e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras, menter_sst;
  spalart_allmaras = neg_spalart_allmaras = e_spalart_allmaras = comp_spalart_allmaras = e_comp_spalart_allmaras = menter_sst = false;

  /*--- Assign turbulence model booleans ---*/

  switch (config->GetKind_Turb_Model()) {
    case TURB_MODEL::NONE:
      SU2_MPI::Error("No turbulence model selected.", CURRENT_FUNCTION);
      break;
    case TURB_MODEL::SA:        spalart_allmaras = true;        break;
    case TURB_MODEL::SA_NEG:    neg_spalart_allmaras = true;    break;
    case TURB_MODEL::SA_E:      e_spalart_allmaras = true;      break;
    case TURB_MODEL::SA_COMP:   comp_spalart_allmaras = true;   break;
    case TURB_MODEL::SA_E_COMP: e_comp_spalart_allmaras = true; break;
    case TURB_MODEL::SST:       menter_sst = true;              break;
    case TURB_MODEL::SST_SUST:  menter_sst = true;              break;
  }

  /*--- If the Menter SST model is used, store the constants of the model and determine the
        free stream values of the turbulent kinetic energy and dissipation rate. ---*/

  const su2double *constants = nullptr;
  su2double kine_Inf = 0.0, omega_Inf = 0.0;

  if (menter_sst) {
    constants = turb_solver->GetConstants();
    kine_Inf  = turb_solver->GetTke_Inf();
    omega_Inf = turb_solver->GetOmega_Inf();
  }

  /*--- Definition of the convective scheme for each equation and mesh level ---*/

  switch (config->GetKind_ConvNumScheme_Turb()) {
    case NO_UPWIND:
      SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_TURB option.", CURRENT_FUNCTION);
      break;
    case SPACE_UPWIND :
      for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        if (spalart_allmaras || neg_spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
          numerics[iMGlevel][TURB_SOL][conv_term] = new CUpwSca_TurbSA<Indices>(nDim, nVar_Turb, config);
        }
        else if (menter_sst) numerics[iMGlevel][TURB_SOL][conv_term] = new CUpwSca_TurbSST<Indices>(nDim, nVar_Turb, config);
      }
      break;
    default:
      SU2_MPI::Error("Invalid convective scheme for the turbulence equations.", CURRENT_FUNCTION);
      break;
  }

  /*--- Definition of the viscous scheme for each equation and mesh level ---*/

  for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
      numerics[iMGlevel][TURB_SOL][visc_term] = new CAvgGrad_TurbSA<Indices>(nDim, nVar_Turb, true, config);
    }
    else if (neg_spalart_allmaras)
      numerics[iMGlevel][TURB_SOL][visc_term] = new CAvgGrad_TurbSA_Neg<Indices>(nDim, nVar_Turb, true, config);
    else if (menter_sst)
      numerics[iMGlevel][TURB_SOL][visc_term] = new CAvgGrad_TurbSST<Indices>(nDim, nVar_Turb, constants, true, config);
  }

  /*--- Definition of the source term integration scheme for each equation and mesh level ---*/

  for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    auto& turb_source_first_term = numerics[iMGlevel][TURB_SOL][source_first_term];

    if (spalart_allmaras)
      turb_source_first_term = new CSourcePieceWise_TurbSA<Indices>(nDim, nVar_Turb, config);
    else if (e_spalart_allmaras)
      turb_source_first_term = new CSourcePieceWise_TurbSA_E<Indices>(nDim, nVar_Turb, config);
    else if (comp_spalart_allmaras)
      turb_source_first_term = new CSourcePieceWise_TurbSA_COMP<Indices>(nDim, nVar_Turb, config);
    else if (e_comp_spalart_allmaras)
      turb_source_first_term = new CSourcePieceWise_TurbSA_E_COMP<Indices>(nDim, nVar_Turb, config);
    else if (neg_spalart_allmaras)
      turb_source_first_term = new CSourcePieceWise_TurbSA_Neg<Indices>(nDim, nVar_Turb, config);
    else if (menter_sst)
      turb_source_first_term = new CSourcePieceWise_TurbSST<Indices>(nDim, nVar_Turb, constants, kine_Inf, omega_Inf,
                                                                     config);

    numerics[iMGlevel][TURB_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Turb, config);
  }

  /*--- Definition of the boundary condition method ---*/

  for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
      numerics[iMGlevel][TURB_SOL][conv_bound_term] = new CUpwSca_TurbSA<Indices>(nDim, nVar_Turb, config);
      numerics[iMGlevel][TURB_SOL][visc_bound_term] = new CAvgGrad_TurbSA<Indices>(nDim, nVar_Turb, false, config);
    }
    else if (neg_spalart_allmaras) {
      numerics[iMGlevel][TURB_SOL][conv_bound_term] = new CUpwSca_TurbSA<Indices>(nDim, nVar_Turb, config);
      numerics[iMGlevel][TURB_SOL][visc_bound_term] = new CAvgGrad_TurbSA_Neg<Indices>(nDim, nVar_Turb, false, config);
    }
    else if (menter_sst) {
      numerics[iMGlevel][TURB_SOL][conv_bound_term] = new CUpwSca_TurbSST<Indices>(nDim, nVar_Turb, config);
      numerics[iMGlevel][TURB_SOL][visc_bound_term] = new CAvgGrad_TurbSST<Indices>(nDim, nVar_Turb, constants, false,
                                                                                    config);
    }
  }
}
/*--- Explicit instantiation of the template above, needed because it is defined in a cpp file, instead of hpp. ---*/
template void CDriver::InstantiateTurbulentNumerics<CEulerVariable::CIndices<unsigned short>>(
    unsigned short, int, const CConfig*, const CSolver*, CNumerics****&) const;

template void CDriver::InstantiateTurbulentNumerics<CIncEulerVariable::CIndices<unsigned short>>(
    unsigned short, int, const CConfig*, const CSolver*, CNumerics****&) const;

template void CDriver::InstantiateTurbulentNumerics<CNEMOEulerVariable::CIndices<unsigned short>>(
    unsigned short, int, const CConfig*, const CSolver*, CNumerics****&) const;

template <class Indices>
void CDriver::InstantiateSpeciesNumerics(unsigned short nVar_Species, int offset, const CConfig *config,
                                         const CSolver* species_solver, CNumerics ****&numerics) const {
  const int conv_term = CONV_TERM + offset;
  const int visc_term = VISC_TERM + offset;

  const int source_first_term = SOURCE_FIRST_TERM + offset;
  const int source_second_term = SOURCE_SECOND_TERM + offset;

  const int conv_bound_term = CONV_BOUND_TERM + offset;
  const int visc_bound_term = VISC_BOUND_TERM + offset;

  /*--- Definition of the convective scheme for each equation and mesh level. Also for boundary conditions. ---*/

  switch (config->GetKind_ConvNumScheme_Species()) {
    case NONE :
      break;
    case SPACE_UPWIND :
      for (auto iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics[iMGlevel][SPECIES_SOL][conv_term] = new CUpwSca_Species<Indices>(nDim, nVar_Species, config);
        numerics[iMGlevel][SPECIES_SOL][conv_bound_term] = new CUpwSca_Species<Indices>(nDim, nVar_Species, config);
      }
      break;
    default :
      SU2_MPI::Error("Invalid convective scheme for the species transport equations. Use SCALAR_UPWIND.", CURRENT_FUNCTION);
      break;
  }

  /*--- Definition of the viscous scheme for each equation and mesh level ---*/

  for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    numerics[iMGlevel][SPECIES_SOL][visc_term] = new CAvgGrad_Species<Indices>(nDim, nVar_Species, true, config);
    numerics[iMGlevel][SPECIES_SOL][visc_bound_term] = new CAvgGrad_Species<Indices>(nDim, nVar_Species, false, config);
  }

  /*--- Definition of the source term integration scheme for each equation and mesh level ---*/

  for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    if (config->GetAxisymmetric() == YES) {
      numerics[iMGlevel][SPECIES_SOL][source_first_term] = new CSourceAxisymmetric_Species<Indices>(nDim, nVar_Species, config);
    }
    else {
      numerics[iMGlevel][SPECIES_SOL][source_first_term] = new CSourceNothing(nDim, nVar_Species, config);
    }
    numerics[iMGlevel][SPECIES_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Species, config);
  }
}
/*--- Explicit instantiation of the template above, needed because it is defined in a cpp file, instead of hpp. ---*/
template void CDriver::InstantiateSpeciesNumerics<CEulerVariable::CIndices<unsigned short>>(
    unsigned short, int, const CConfig*, const CSolver*, CNumerics****&) const;

template void CDriver::InstantiateSpeciesNumerics<CIncEulerVariable::CIndices<unsigned short>>(
    unsigned short, int, const CConfig*, const CSolver*, CNumerics****&) const;

template void CDriver::InstantiateSpeciesNumerics<CNEMOEulerVariable::CIndices<unsigned short>>(
    unsigned short, int, const CConfig*, const CSolver*, CNumerics****&) const;

void CDriver::Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics) const {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Numerics Preprocessing ( Zone " << config->GetiZone() <<" ) -------------------" << endl;

  unsigned short iMGlevel, iSol,

  nVar_Template         = 0,
  nVar_Flow             = 0,
  nVar_NEMO             = 0,
  nPrimVar_NEMO         = 0,
  nPrimVarGrad_NEMO     = 0,
  nVar_Trans            = 0,
  nVar_Turb             = 0,
  nVar_Species          = 0,
  nVar_Adj_Flow         = 0,
  nVar_Adj_Turb         = 0,
  nVar_FEM              = 0,
  nVar_Rad              = 0,
  nVar_Heat             = 0;

  numerics = new CNumerics***[config->GetnMGLevels()+1] ();

  bool compressible = false;
  bool incompressible = false;
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR) || (config->GetKind_FluidModel() == IDEAL_GAS);
  bool roe_low_dissipation = (config->GetKind_RoeLowDiss() != NO_ROELOWDISS);

  /*--- Initialize some useful booleans ---*/
  bool euler, ns, NEMO_euler, NEMO_ns, turbulent, adj_euler, adj_ns, adj_turb, fem_euler, fem_ns;
  bool fem, heat, transition, template_solver;

  euler = ns = NEMO_euler = NEMO_ns = turbulent = adj_euler = adj_ns = adj_turb = fem_euler = fem_ns = false;
  fem = heat = transition = template_solver = false;
  bool species = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case MAIN_SOLVER::TEMPLATE_SOLVER:
      template_solver = true; break;

    case MAIN_SOLVER::EULER:
    case MAIN_SOLVER::DISC_ADJ_EULER:
      euler = compressible = true; break;

    case MAIN_SOLVER::NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
      ns = compressible = true;
      species = (config->GetKind_Species_Model() != SPECIES_MODEL::NONE); break;

    case MAIN_SOLVER::NEMO_EULER:
      NEMO_euler = compressible = true; break;

    case MAIN_SOLVER::NEMO_NAVIER_STOKES:
      NEMO_ns = compressible = true; break;

    case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::DISC_ADJ_RANS:
      ns = compressible = turbulent = true;
      transition = (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM);
      species = config->GetKind_Species_Model() != SPECIES_MODEL::NONE; break;

    case MAIN_SOLVER::INC_EULER:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER:
      euler = incompressible = true; break;

    case MAIN_SOLVER::INC_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:
      ns = incompressible = true;
      heat = config->GetWeakly_Coupled_Heat();
      species = (config->GetKind_Species_Model() != SPECIES_MODEL::NONE); break;

    case MAIN_SOLVER::INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      ns = incompressible = turbulent = true;
      heat = config->GetWeakly_Coupled_Heat();
      transition = (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM);
      species = (config->GetKind_Species_Model() != SPECIES_MODEL::NONE); break;

    case MAIN_SOLVER::FEM_EULER:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:
      fem_euler = compressible = true; break;

    case MAIN_SOLVER::FEM_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_FEM_NS:
      fem_ns = compressible = true; break;

    case MAIN_SOLVER::FEM_RANS:
    case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
      fem_ns = compressible = true; break;

    case MAIN_SOLVER::FEM_LES:
      fem_ns = compressible = true; break;

    case MAIN_SOLVER::HEAT_EQUATION:
    case MAIN_SOLVER::DISC_ADJ_HEAT:
      heat = true; break;

    case MAIN_SOLVER::FEM_ELASTICITY:
    case MAIN_SOLVER::DISC_ADJ_FEM:
      fem = true; break;

    case MAIN_SOLVER::ADJ_EULER:
      adj_euler = euler = compressible = true; break;

    case MAIN_SOLVER::ADJ_NAVIER_STOKES:
      adj_ns = ns = compressible = true;
      turbulent = (config->GetKind_Turb_Model() != TURB_MODEL::NONE); break;

    case MAIN_SOLVER::ADJ_RANS:
      adj_ns = ns = compressible = turbulent = true;
      adj_turb = !config->GetFrozen_Visc_Cont(); break;

    default:
      break;

  }

  /*--- Number of variables for the template ---*/

  if (template_solver) nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();

  /*--- Number of variables for direct problem ---*/

  if (euler)        nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)           nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (NEMO_euler)   nVar_NEMO = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (NEMO_ns)      nVar_NEMO = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)    nVar_Turb = solver[MESH_0][TURB_SOL]->GetnVar();
  if (transition)   nVar_Trans = solver[MESH_0][TRANS_SOL]->GetnVar();
  if (species)      nVar_Species = solver[MESH_0][SPECIES_SOL]->GetnVar();

  if (fem_euler)    nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (fem_ns)       nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();

  if (fem)          nVar_FEM = solver[MESH_0][FEA_SOL]->GetnVar();
  if (heat)         nVar_Heat = solver[MESH_0][HEAT_SOL]->GetnVar();

  if (config->AddRadiation()) nVar_Rad = solver[MESH_0][RAD_SOL]->GetnVar();

  /*--- Number of variables for adjoint problem ---*/

  if (adj_euler)    nVar_Adj_Flow = solver[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)       nVar_Adj_Flow = solver[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)     nVar_Adj_Turb = solver[MESH_0][ADJTURB_SOL]->GetnVar();

  /*--- Additional Variables required for NEMO solver ---*/

  if (NEMO_euler || NEMO_ns) nPrimVar_NEMO     = solver[MESH_0][FLOW_SOL]->GetnPrimVar();
  if (NEMO_euler || NEMO_ns) nPrimVarGrad_NEMO = solver[MESH_0][FLOW_SOL]->GetnPrimVarGrad();

  /*--- Definition of the Class for the numerical method:
    numerics_container[INSTANCE_LEVEL][MESH_LEVEL][EQUATION][EQ_TERM] ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    numerics[iMGlevel] = new CNumerics** [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      numerics[iMGlevel][iSol] = new CNumerics* [MAX_TERMS*omp_get_max_threads()]();
  }

  /*--- Instantiate one numerics object per thread for each required term. ---*/

  for (int thread = 0; thread < omp_get_max_threads(); ++thread)
  {
  const int offset = thread * MAX_TERMS;

  const int conv_term = CONV_TERM + offset;
  const int visc_term = VISC_TERM + offset;

  const int source_first_term = SOURCE_FIRST_TERM + offset;
  const int source_second_term = SOURCE_SECOND_TERM + offset;

  const int conv_bound_term = CONV_BOUND_TERM + offset;
  const int visc_bound_term = VISC_BOUND_TERM + offset;

  const int fea_term = FEA_TERM + offset;

  /*--- Solver definition for the template problem ---*/
  if (template_solver) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][TEMPLATE_SOL][conv_term] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default:
        SU2_MPI::Error("Convective scheme not implemented (template_solver).", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][TEMPLATE_SOL][visc_term] = new CViscous_Template(nDim, nVar_Template, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][TEMPLATE_SOL][source_first_term] = new CSource_Template(nDim, nVar_Template, config);

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TEMPLATE_SOL][conv_bound_term] = new CConvective_Template(nDim, nVar_Template, config);
    }

  }

  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_FLOW option.", CURRENT_FUNCTION);
        break;

      case SPACE_CENTERED :
        if (compressible) {
          /*--- "conv_term" is not instantiated as all compressible centered schemes are vectorized. ---*/

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);

        }
        if (incompressible) {
          /*--- Incompressible flow, use preconditioning method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentLaxInc_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentJSTInc_Flow(nDim, nVar_Flow, config); break;
            default:
              SU2_MPI::Error("Invalid centered scheme or not implemented.\n Currently, only JST and LAX-FRIEDRICH are available for incompressible flows.", CURRENT_FUNCTION);
              break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_term] = new CCentLaxInc_Flow(nDim, nVar_Flow, config);

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);

        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE:
              if (ideal_gas) {

                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);
                }
              } else {

                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                }
              }
              break;

            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;

            case AUSMPLUSUP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
              }
              break;

            case AUSMPLUSUP2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
              }
              break;

            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;

            case L2ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
              }
              break;
            case LMROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
              }
              break;

            case SLAU:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, false);
              }
              break;

            case SLAU2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, false);
              }
              break;

            case HLLC:
              if (ideal_gas) {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              else {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              break;

            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;

            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;

            default:
              SU2_MPI::Error("Invalid upwind scheme or not implemented.", CURRENT_FUNCTION);
              break;
          }

        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case FDS:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
              }
              break;
            default:
              SU2_MPI::Error("Invalid upwind scheme or not implemented.\n Currently, only FDS is available for incompressible flows.", CURRENT_FUNCTION);
              break;
          }
        }
        break;

      default:
        SU2_MPI::Error("Invalid convective scheme for the Euler / Navier-Stokes equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible) {
      if (ideal_gas) {

        /*--- Compressible flow Ideal gas ---*/
        numerics[MESH_0][FLOW_SOL][visc_term] = new CAvgGrad_Flow(nDim, nVar_Flow, true, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_term] = new CAvgGrad_Flow(nDim, nVar_Flow, false, config);

        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CAvgGrad_Flow(nDim, nVar_Flow, false, config);

      } else {

        /*--- Compressible flow Real gas ---*/
        numerics[MESH_0][FLOW_SOL][visc_term] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, true, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_term] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, false, config);

        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, false, config);

      }
    }
    if (incompressible) {
      /*--- Incompressible flow, use preconditioning method ---*/
      numerics[MESH_0][FLOW_SOL][visc_term] = new CAvgGradInc_Flow(nDim, nVar_Flow, true, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics[iMGlevel][FLOW_SOL][visc_term] = new CAvgGradInc_Flow(nDim, nVar_Flow, false, config);

      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CAvgGradInc_Flow(nDim, nVar_Flow, false, config);
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      if (config->GetBody_Force() == YES) {
        if (incompressible)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncBodyForce(nDim, nVar_Flow, config);
        else
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceBodyForce(nDim, nVar_Flow, config);
      }
      else if (incompressible && (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncStreamwise_Periodic(nDim, nVar_Flow, config);
      }
      else if (incompressible && (config->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ)) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceBoussinesq(nDim, nVar_Flow, config);
      }
      else if (config->GetRotating_Frame() == YES) {
        if (incompressible)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncRotatingFrame_Flow(nDim, nVar_Flow, config);
        else
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
      }
      else if (config->GetAxisymmetric() == YES) {
        if (incompressible)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncAxisymmetric_Flow(nDim, nVar_Flow, config);
        else if (ideal_gas)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceAxisymmetric_Flow(nDim, nVar_Flow, config);
        else
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceGeneralAxisymmetric_Flow(nDim, nVar_Flow, config);
      }
      else if (config->GetGravityForce() == YES) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceGravity(nDim, nVar_Flow, config);
      }
      else if (config->GetWind_Gust() == YES) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceWindGust(nDim, nVar_Flow, config);
      }
      else {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceNothing(nDim, nVar_Flow, config);
      }

      /*--- At the moment it is necessary to have the RHT equation in order to have a volumetric heat source. ---*/
      if (config->AddRadiation())
        numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceRadiation(nDim, nVar_Flow, config);
      else if ((incompressible && (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)) &&
               (config->GetEnergy_Equation() && !config->GetStreamwise_Periodic_Temperature()))
        numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceIncStreamwisePeriodic_Outlet(nDim, nVar_Flow, config);
      else
        numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Flow, config);
    }

  }

   /*--- Solver definition for the Potential, Euler, Navier-Stokes NEMO problems ---*/

  if (NEMO_euler || NEMO_ns) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_FLOW option.", CURRENT_FUNCTION);
        break;

      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentLax_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config); break;
            default:
            SU2_MPI::Error("Invalid centered scheme or not implemented.", CURRENT_FUNCTION);
            break;
          }

          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_term] = new CCentLax_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwRoe_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
              }
              break;

            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSM_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSM_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
              }
              break;

            case AUSMPLUSUP2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP2_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP2_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
              }
              break;

            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
              }
              break;

            case AUSMPWPLUS:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPWplus_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPWplus_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
              }
              break;

            default:
              SU2_MPI::Error("Invalid upwind scheme or not implemented.", CURRENT_FUNCTION);
              break;
          }

        }
        break;

      default:
        SU2_MPI::Error("Invalid convective scheme for the NEMO Euler / Navier-Stokes equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible) {

      numerics[MESH_0][FLOW_SOL][visc_term] = new CAvgGradCorrected_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics[iMGlevel][FLOW_SOL][visc_term] = new CAvgGrad_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);

      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CAvgGrad_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSource_NEMO(nDim, nVar_NEMO, nPrimVar_NEMO, nPrimVarGrad_NEMO, config);
      numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceNothing(nDim, nVar_NEMO, config);
    }
  }

  /*--- Riemann solver definition for the Euler, Navier-Stokes problems for the FEM discretization. ---*/
  if ((fem_euler) || (fem_ns)) {

    switch (config->GetRiemann_Solver_FEM()) {
      case ROE:
      case LAX_FRIEDRICH:
        /* Hard coded optimized implementation is used in the DG solver. No need to allocate the
           corresponding entry in numerics. */
        break;

      case AUSM:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
        }
        break;

      case TURKEL:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
        }
        break;

      case HLLC:
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
            numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
          }
        break;

      case MSW:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
        }
        break;

      case CUSP:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
        }
        break;

      default:
        SU2_MPI::Error("Riemann solver not implemented.", CURRENT_FUNCTION);
        break;
    }

  }

  /*--- Solver definition for the turbulent model problem ---*/

  if (turbulent) {
    if (incompressible)
      InstantiateTurbulentNumerics<CIncEulerVariable::CIndices<unsigned short> >(nVar_Turb, offset, config,
                                                                                 solver[MESH_0][TURB_SOL], numerics);
    else if (NEMO_ns)
      InstantiateTurbulentNumerics<CNEMOEulerVariable::CIndices<unsigned short> >(nVar_Turb, offset, config,
                                                                                  solver[MESH_0][TURB_SOL], numerics);
    else
      InstantiateTurbulentNumerics<CEulerVariable::CIndices<unsigned short> >(nVar_Turb, offset, config,
                                                                              solver[MESH_0][TURB_SOL], numerics);
  }

  /*--- Solver definition for the transition model problem ---*/
  if (transition) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NO_UPWIND:
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_TURB option.", CURRENT_FUNCTION);
        break;
      case SPACE_UPWIND:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][TRANS_SOL][conv_term] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
        }
        break;
      default:
        SU2_MPI::Error("Invalid convective scheme for the transition equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TRANS_SOL][visc_term] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TRANS_SOL][source_first_term] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
      numerics[iMGlevel][TRANS_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Trans, config);
    }

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TRANS_SOL][conv_bound_term] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
    }
  }

  /*--- Solver definition for the species transport problem ---*/

  if (species) {
    if (incompressible)
      InstantiateSpeciesNumerics<CIncEulerVariable::CIndices<unsigned short> >(nVar_Species, offset, config,
                                                                               solver[MESH_0][SPECIES_SOL], numerics);
    else if (compressible)
      InstantiateSpeciesNumerics<CEulerVariable::CIndices<unsigned short> >(nVar_Species, offset, config,
                                                                            solver[MESH_0][SPECIES_SOL], numerics);
    else
      SU2_MPI::Error("Species transport only available for standard compressible and incompressible flow.", CURRENT_FUNCTION);
  }

  /*--- Solver definition of the finite volume heat solver  ---*/
  if (heat) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      numerics[iMGlevel][HEAT_SOL][visc_term] = new CAvgGrad_Heat(nDim, nVar_Heat, config, true);
      numerics[iMGlevel][HEAT_SOL][visc_bound_term] = new CAvgGrad_Heat(nDim, nVar_Heat, config, false);

      switch (config->GetKind_ConvNumScheme_Heat()) {

        case SPACE_UPWIND :
          numerics[iMGlevel][HEAT_SOL][conv_term] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          numerics[iMGlevel][HEAT_SOL][conv_bound_term] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          break;

        case SPACE_CENTERED :
          numerics[iMGlevel][HEAT_SOL][conv_term] = new CCentSca_Heat(nDim, nVar_Heat, config);
          numerics[iMGlevel][HEAT_SOL][conv_bound_term] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          break;

        default:
          SU2_MPI::Error("Invalid convective scheme for the heat transfer equations.", CURRENT_FUNCTION);
          break;
      }
    }
  }

  /*--- Solver definition for the radiation model problem ---*/

  if (config->AddRadiation()) {
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics[MESH_0][RAD_SOL][VISC_TERM] = new CAvgGradCorrected_P1(nDim, nVar_Rad, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics[MESH_0][RAD_SOL][SOURCE_FIRST_TERM] = new CSourceP1(nDim, nVar_Rad, config);

    /*--- Definition of the boundary condition method ---*/
    numerics[MESH_0][RAD_SOL][VISC_BOUND_TERM] = new CAvgGradCorrected_P1(nDim, nVar_Rad, config);
  }

  /*--- Solver definition for the flow adjoint problem ---*/

  if (adj_euler || adj_ns) {

    if (incompressible)
      SU2_MPI::Error("Convective schemes not implemented for incompressible continuous adjoint.", CURRENT_FUNCTION);

    /*--- Definition of the convective scheme for each equation and mesh level ---*/

    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE:
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_ADJFLOW option.", CURRENT_FUNCTION);
        break;

      case SPACE_CENTERED :

        if (compressible) {

          /*--- Compressible flow ---*/

          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : numerics[MESH_0][ADJFLOW_SOL][conv_term] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics[MESH_0][ADJFLOW_SOL][conv_term] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default:
              SU2_MPI::Error("Centered scheme not implemented.", CURRENT_FUNCTION);
              break;
          }

          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][ADJFLOW_SOL][conv_term] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);

          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][ADJFLOW_SOL][conv_bound_term] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);

        }
        break;

      case SPACE_UPWIND :

        if (compressible) {

          /*--- Compressible flow ---*/

          switch (config->GetKind_Upwind_AdjFlow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][ADJFLOW_SOL][conv_term] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics[iMGlevel][ADJFLOW_SOL][conv_bound_term] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default:
              SU2_MPI::Error("Upwind scheme not implemented.", CURRENT_FUNCTION);
              break;
          }
        }
        break;

      default:
        SU2_MPI::Error("Invalid convective scheme for the continuous adjoint Euler / Navier-Stokes equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/

    if (compressible) {

      /*--- Compressible flow ---*/

      numerics[MESH_0][ADJFLOW_SOL][visc_term] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics[MESH_0][ADJFLOW_SOL][visc_bound_term] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);

      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics[iMGlevel][ADJFLOW_SOL][visc_term] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics[iMGlevel][ADJFLOW_SOL][visc_bound_term] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      }

    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/

      if (compressible) {

        if (adj_ns) {

          numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceViscous_AdjFlow(nDim, nVar_Adj_Flow, config);

          if (config->GetRotating_Frame() == YES)
            numerics[iMGlevel][ADJFLOW_SOL][source_second_term] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics[iMGlevel][ADJFLOW_SOL][source_second_term] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);

        }

        else {

          if (config->GetRotating_Frame() == YES)
            numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else if (config->GetAxisymmetric() == YES)
            numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceAxisymmetric_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceNothing(nDim, nVar_Adj_Flow, config);

          numerics[iMGlevel][ADJFLOW_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Adj_Flow, config);

        }

      }

    }

  }

  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {

    if (config->GetKind_Turb_Model() != TURB_MODEL::SA)
      SU2_MPI::Error("Only the SA turbulence model can be used with the continuous adjoint solver.", CURRENT_FUNCTION);

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {
      case NO_CONVECTIVE:
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_ADJTURB option.", CURRENT_FUNCTION);
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][ADJTURB_SOL][conv_term] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
        break;
      default:
        SU2_MPI::Error("Convective scheme not implemented (adjoint turbulence).", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][ADJTURB_SOL][visc_term] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][ADJTURB_SOL][source_first_term] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
      numerics[iMGlevel][ADJTURB_SOL][source_second_term] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
    }

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][ADJTURB_SOL][conv_bound_term] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);

  }

  /*--- Numerics definition for FEM-like problems. ---*/

  if (fem) {
    /*--- Initialize the container for FEA_TERM. This will be the only one for most of the cases. ---*/
    switch (config->GetGeometricConditions()) {
      case STRUCT_DEFORMATION::SMALL:
        switch (config->GetMaterialModel()) {
          case STRUCT_MODEL::LINEAR_ELASTIC:
            numerics[MESH_0][FEA_SOL][fea_term] = new CFEALinearElasticity(nDim, nVar_FEM, config);
            break;
          default:
            SU2_MPI::Error("Material model does not correspond to geometric conditions.", CURRENT_FUNCTION);
            break;
        }
        break;
      case STRUCT_DEFORMATION::LARGE:
        switch (config->GetMaterialModel()) {
          case STRUCT_MODEL::LINEAR_ELASTIC:
            SU2_MPI::Error("Material model does not correspond to geometric conditions.", CURRENT_FUNCTION);
            break;
          case STRUCT_MODEL::NEO_HOOKEAN:
            if (config->GetMaterialCompressibility() == STRUCT_COMPRESS::COMPRESSIBLE) {
              numerics[MESH_0][FEA_SOL][fea_term] = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config);
            } else {
              SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            }
            break;
          case STRUCT_MODEL::KNOWLES:
            if (config->GetMaterialCompressibility() == STRUCT_COMPRESS::NEARLY_INCOMP) {
              numerics[MESH_0][FEA_SOL][fea_term] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config);
            } else {
              SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            }
            break;
          case STRUCT_MODEL::IDEAL_DE:
            if (config->GetMaterialCompressibility() == STRUCT_COMPRESS::NEARLY_INCOMP) {
              numerics[MESH_0][FEA_SOL][fea_term] = new CFEM_IdealDE(nDim, nVar_FEM, config);
            } else {
              SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            }
            break;
        }
        break;
    }

    /*--- The following definitions only make sense if we have a non-linear solution. ---*/
    if (config->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE) {

      /*--- This allocates a container for electromechanical effects. ---*/

      bool de_effects = config->GetDE_Effects();
      if (de_effects)
        numerics[MESH_0][FEA_SOL][DE_TERM+offset] = new CFEM_DielectricElastomer(nDim, nVar_FEM, config);

      ifstream properties_file;

      string filename = config->GetFEA_FileName();
      if (nZone > 1)
        filename = config->GetMultizone_FileName(filename, iZone, ".dat");

      properties_file.open(filename.data(), ios::in);

      /*--- In case there is a properties file, containers are allocated for a number of material models. ---*/

      if (!(properties_file.fail())) {
        numerics[MESH_0][FEA_SOL][MAT_NHCOMP+offset]  = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config);
        numerics[MESH_0][FEA_SOL][MAT_IDEALDE+offset] = new CFEM_IdealDE(nDim, nVar_FEM, config);
        numerics[MESH_0][FEA_SOL][MAT_KNOWLES+offset] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config);
      }
    }
  }

  /*--- Instantiate the numerics for the mesh solver. ---*/
  if (config->GetDeform_Mesh())
    numerics[MESH_0][MESH_SOL][fea_term] = new CFEAMeshElasticity(nDim, nDim, geometry[MESH_0]->GetnElem(), config);

  } // end "per-thread" allocation loop

}

void CDriver::Numerics_Postprocessing(CNumerics *****numerics, CSolver***, CGeometry**,
                                      CConfig *config, unsigned short val_iInst) {

  for (unsigned short iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++) {

      for (unsigned int iTerm = 0; iTerm < MAX_TERMS*omp_get_max_threads(); iTerm++) {

        delete numerics[val_iInst][iMGlevel][iSol][iTerm];
      }
      delete [] numerics[val_iInst][iMGlevel][iSol];
    }
    delete[] numerics[val_iInst][iMGlevel];
  }
  delete[] numerics[val_iInst];

}

void CDriver::Iteration_Preprocessing(CConfig* config, CIteration *&iteration) const {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Iteration Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  iteration = CIterationFactory::CreateIteration(config->GetKind_Solver(), config);

}

void CDriver::DynamicMesh_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CIteration* iteration,
                                        CVolumetricMovement *&grid_movement, CSurfaceMovement *&surface_movement) const{

  /*--- Instantiate the geometry movement classes for the solution of unsteady
   flows on dynamic meshes, including rigid mesh transformations, dynamically
   deforming meshes, and preprocessing of harmonic balance. ---*/

  if (!fem_solver && (config->GetGrid_Movement() || (config->GetDirectDiff() == D_DESIGN))) {
    if (rank == MASTER_NODE)
      cout << "Setting dynamic mesh structure for zone "<< iZone + 1<<"." << endl;
    grid_movement = new CVolumetricMovement(geometry[MESH_0], config);

    surface_movement = new CSurfaceMovement();
    surface_movement->CopyBoundary(geometry[MESH_0], config);
    if (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE){
      if (rank == MASTER_NODE) cout << endl <<  "Instance "<< iInst + 1 <<":" << endl;
      iteration->SetGrid_Movement(geometry, surface_movement, grid_movement,  solver, config, 0, iInst);
    }
  }

  if (config->GetDirectDiff() == D_DESIGN) {
    if (rank == MASTER_NODE)
      cout << "Setting surface/volume derivatives." << endl;

    /*--- Set the surface derivatives, i.e. the derivative of the surface mesh nodes with respect to the design variables ---*/

    surface_movement->SetSurface_Derivative(geometry[MESH_0],config);

    /*--- Call the volume deformation routine with derivative mode enabled.
       This computes the derivative of the volume mesh with respect to the surface nodes ---*/

    grid_movement->SetVolume_Deformation(geometry[MESH_0],config, true, true);

    /*--- Update the multi-grid structure to propagate the derivative information to the coarser levels ---*/

    CGeometry::UpdateGeometry(geometry,config);

  }

}

void CDriver::Interface_Preprocessing(CConfig **config, CSolver***** solver, CGeometry**** geometry,
                                      unsigned short** interface_types, CInterface ***interface,
                                      vector<vector<unique_ptr<CInterpolator> > >& interpolation) {

  /*--- Setup interpolation and transfer for all possible donor/target pairs. ---*/

  for (auto target = 0u; target < nZone; target++) {

    for (auto donor = 0u; donor < nZone; donor++) {

      /*--- Aliases to make code less verbose. ---*/
      auto& interface_type = interface_types[donor][target];

      if (donor == target) {
        interface_type = ZONES_ARE_EQUAL;
        continue;
      }
      interface_type = NO_TRANSFER;

      /*--- If there is a common interface setup the interpolation and transfer. ---*/

      if (!CInterpolator::CheckZonesInterface(config[donor], config[target])) {
        interface_type = NO_COMMON_INTERFACE;
      }
      else {
        /*--- Begin the creation of the communication pattern among zones. ---*/

        if (rank == MASTER_NODE) cout << "From zone " << donor << " to zone " << target << ":" << endl;

        /*--- Setup the interpolation. ---*/

        interpolation[donor][target] = unique_ptr<CInterpolator>(CInterpolatorFactory::CreateInterpolator(
                                       geometry, config, interpolation[target][donor].get(), donor, target));

        /*--- The type of variables transferred depends on the donor/target physics. ---*/

        const bool heat_target = config[target]->GetHeatProblem();
        const bool fluid_target = config[target]->GetFluidProblem();
        const bool structural_target = config[target]->GetStructuralProblem();

        const bool heat_donor = config[donor]->GetHeatProblem();
        const bool fluid_donor = config[donor]->GetFluidProblem();
        const bool structural_donor = config[donor]->GetStructuralProblem();

        /*--- Initialize the appropriate transfer strategy. ---*/

        if (rank == MASTER_NODE) cout << " Transferring ";

        if (fluid_donor && structural_target) {
          interface_type = FLOW_TRACTION;
          auto nConst = 2;
          bool conservative = config[target]->GetConservativeInterpolation();
          if(!config[ZONE_0]->GetDiscrete_Adjoint()) {
            interface[donor][target] = new CFlowTractionInterface(nDim, nConst, config[donor], conservative);
          } else {
            interface[donor][target] = new CDiscAdjFlowTractionInterface(nDim, nConst, config[donor], conservative);
          }
          if (rank == MASTER_NODE) cout << "fluid " << (conservative? "forces." : "tractions.") << endl;
        }
        else if (structural_donor && (fluid_target || heat_target)) {
          if (solver_container[target][INST_0][MESH_0][MESH_SOL] == nullptr) {
            SU2_MPI::Error("Mesh deformation was not correctly specified for the fluid/heat zone.\n"
                           "Use DEFORM_MESH=YES, and setup MARKER_DEFORM_MESH=(...)", CURRENT_FUNCTION);
          }
          interface_type = BOUNDARY_DISPLACEMENTS;
          if (!config[donor]->GetTime_Domain()) interface[donor][target] = new CDisplacementsInterface(nDim, 0);
          else interface[donor][target] = new CDisplacementsInterface(2*nDim, 0);
          if (rank == MASTER_NODE) cout << "boundary displacements from the structural solver." << endl;
        }
        else if (fluid_donor && fluid_target) {
          interface_type = SLIDING_INTERFACE;
          auto nVar = solver[donor][INST_0][MESH_0][FLOW_SOL]->GetnPrimVar();
          interface[donor][target] = new CSlidingInterface(nVar, 0);
          if (rank == MASTER_NODE) cout << "sliding interface." << endl;
        }
        else if (heat_donor || heat_target) {
          if (heat_donor && heat_target)
            SU2_MPI::Error("Conjugate heat transfer between solids is not implemented.", CURRENT_FUNCTION);

          const auto fluidZone = heat_target? donor : target;

          if (config[fluidZone]->GetEnergy_Equation() || (config[fluidZone]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE))
            interface_type = heat_target? CONJUGATE_HEAT_FS : CONJUGATE_HEAT_SF;
          else if (config[fluidZone]->GetWeakly_Coupled_Heat())
            interface_type = heat_target? CONJUGATE_HEAT_WEAKLY_FS : CONJUGATE_HEAT_WEAKLY_SF;
          else
            interface_type = NO_TRANSFER;

          if (interface_type != NO_TRANSFER) {
            auto nVar = 4;
            interface[donor][target] = new CConjugateHeatInterface(nVar, 0);
            if (rank == MASTER_NODE) cout << "conjugate heat variables." << endl;
          }
          else {
            if (rank == MASTER_NODE) cout << "NO heat variables." << endl;
          }
        }
        else {
          if (solver[donor][INST_0][MESH_0][FLOW_SOL] == nullptr)
            SU2_MPI::Error("Could not determine the number of variables for transfer.", CURRENT_FUNCTION);

          auto nVar = solver[donor][INST_0][MESH_0][FLOW_SOL]->GetnVar();
          interface_type = CONSERVATIVE_VARIABLES;
          interface[donor][target] = new CConservativeVarsInterface(nVar, 0);
          if (rank == MASTER_NODE) cout << "generic conservative variables." << endl;
        }
      }

      /*--- Mixing plane for turbo machinery applications. ---*/

      if (config[donor]->GetBoolMixingPlaneInterface()) {
        interface_type = MIXING_PLANE;
        auto nVar = solver[donor][INST_0][MESH_0][FLOW_SOL]->GetnVar();
        interface[donor][target] = new CMixingPlaneInterface(nVar, 0);
        if (rank == MASTER_NODE) {
          cout << "Set mixing-plane interface from donor zone "
               << donor << " to target zone " << target << "." << endl;
        }
      }

    }

  }

}

void CDriver::StaticMesh_Preprocessing(const CConfig *config, CGeometry** geometry){

  unsigned short iMGlevel, iMGfine;
  unsigned short Kind_Grid_Movement;

  unsigned short iZone = config->GetiZone();

  Kind_Grid_Movement = config->GetKind_GridMovement();

  if (!fem_solver) {

    switch (Kind_Grid_Movement) {

      case ROTATING_FRAME:

        /*--- Steadily rotating frame: set the grid velocities just once
         before the first iteration flow solver. ---*/

        if (rank == MASTER_NODE) {
          cout << endl << " Setting rotating frame grid velocities";
          cout << " for zone " << iZone << "." << endl;
        }

        /*--- Set the grid velocities on all multigrid levels for a steadily
           rotating reference frame. ---*/

        for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++){
          geometry[iMGlevel]->SetRotationalVelocity(config, true);
          geometry[iMGlevel]->SetShroudVelocity(config);
        }

        break;

      case STEADY_TRANSLATION:

        /*--- Set the translational velocity and hold the grid fixed during
         the calculation (similar to rotating frame, but there is no extra
         source term for translation). ---*/

        if (rank == MASTER_NODE)
          cout << endl << " Setting translational grid velocities." << endl;

        /*--- Set the translational velocity on all grid levels. ---*/

        for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++)
          geometry_container[iZone][INST_0][iMGlevel]->SetTranslationalVelocity(config, true);

        break;

      default:
        break;
    }

    if (config->GetnMarker_Moving() > 0) {

      /*--- Fixed wall velocities: set the grid velocities only one time
       before the first iteration flow solver. ---*/
      if (rank == MASTER_NODE)
        cout << endl << " Setting the moving wall velocities." << endl;

      geometry[MESH_0]->SetWallVelocity(config, true);

      /*--- Update the grid velocities on the coarser multigrid levels after
        setting the moving wall velocities for the finest mesh. ---*/
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++){
        iMGfine = iMGlevel-1;
        geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine]);
      }
    }
  } else {

    /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
         define all virtual functions in the base class CGeometry. ---*/
    CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry[MESH_0]);

    /*--- Initialize the static mesh movement, if necessary. ---*/
    const unsigned short Kind_Grid_Movement = config->GetKind_GridMovement();
    const bool initStaticMovement = (config->GetGrid_Movement() &&
                                     (Kind_Grid_Movement == MOVING_WALL    ||
                                      Kind_Grid_Movement == ROTATING_FRAME ||
                                      Kind_Grid_Movement == STEADY_TRANSLATION));

    if(initStaticMovement){
      if (rank == MASTER_NODE) cout << "Initialize Static Mesh Movement" << endl;
      DGMesh->InitStaticMeshMovement(config, Kind_Grid_Movement, iZone);
    }
  }

}

void CDriver::Output_Preprocessing(CConfig **config, CConfig *driver_config, COutput **&output, COutput *&driver_output){

  /*--- Definition of the output class (one for each zone). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    if (rank == MASTER_NODE)
      cout << endl <<"-------------------- Output Preprocessing ( Zone " << iZone <<" ) --------------------" << endl;

    MAIN_SOLVER kindSolver = config[iZone]->GetKind_Solver();

    output[iZone] = COutputFactory::CreateOutput(kindSolver, config[iZone], nDim);

    /*--- If dry-run is used, do not open/overwrite history file. ---*/
    output[iZone]->PreprocessHistoryOutput(config[iZone], !dry_run);

    output[iZone]->PreprocessVolumeOutput(config[iZone]);

  }

  if (driver_config->GetMultizone_Problem()){
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Output Preprocessing ( Multizone ) ------------------" << endl;

    driver_output = COutputFactory::CreateMultizoneOutput(driver_config, config, nDim);

    driver_output->PreprocessMultizoneHistoryOutput(output, config, driver_config, !dry_run);
  }

  /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetTime_Domain() && config_container[ZONE_0]->GetRestart())
    TimeIter = config_container[ZONE_0]->GetRestart_Iter();

}


void CDriver::Turbomachinery_Preprocessing(CConfig** config, CGeometry**** geometry, CSolver***** solver,
                                           CInterface*** interface){

  unsigned short donorZone,targetZone, nMarkerInt, iMarkerInt;
  unsigned short nSpanMax = 0;
  bool restart   = (config[ZONE_0]->GetRestart() || config[ZONE_0]->GetRestart_Flow());
  mixingplane = config[ZONE_0]->GetBoolMixingPlaneInterface();
  bool discrete_adjoint = config[ZONE_0]->GetDiscrete_Adjoint();
  su2double areaIn, areaOut, nBlades, flowAngleIn, flowAngleOut;

  /*--- Create turbovertex structure ---*/
  if (rank == MASTER_NODE) cout<<endl<<"Initialize Turbo Vertex Structure." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config[iZone]->GetBoolTurbomachinery()){
      geometry[iZone][INST_0][MESH_0]->ComputeNSpan(config[iZone], iZone, INFLOW, true);
      geometry[iZone][INST_0][MESH_0]->ComputeNSpan(config[iZone], iZone, OUTFLOW, true);
      if (rank == MASTER_NODE) cout <<"Number of span-wise sections in Zone "<< iZone<<": "<< config[iZone]->GetnSpanWiseSections() <<"."<< endl;
      if (config[iZone]->GetnSpanWiseSections() > nSpanMax){
        nSpanMax = config[iZone]->GetnSpanWiseSections();
      }

      config[ZONE_0]->SetnSpan_iZones(config[iZone]->GetnSpanWiseSections(), iZone);

      geometry[iZone][INST_0][MESH_0]->SetTurboVertex(config[iZone], iZone, INFLOW, true);
      geometry[iZone][INST_0][MESH_0]->SetTurboVertex(config[iZone], iZone, OUTFLOW, true);
    }
  }

  /*--- Set maximum number of Span among all zones ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config[iZone]->GetBoolTurbomachinery()){
      config[iZone]->SetnSpanMaxAllZones(nSpanMax);
    }
  }
  if (rank == MASTER_NODE) cout<<"Max number of span-wise sections among all zones: "<< nSpanMax<<"."<< endl;


  if (rank == MASTER_NODE) cout<<"Initialize solver containers for average and performance quantities." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->InitTurboContainers(geometry[iZone][INST_0][MESH_0],config[iZone]);
  }

//TODO(turbo) make it general for turbo HB
  if (rank == MASTER_NODE) cout<<"Compute inflow and outflow average geometric quantities." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    geometry[iZone][INST_0][MESH_0]->SetAvgTurboValue(config[iZone], iZone, INFLOW, true);
    geometry[iZone][INST_0][MESH_0]->SetAvgTurboValue(config[iZone],iZone, OUTFLOW, true);
    geometry[iZone][INST_0][MESH_0]->GatherInOutAverageValues(config[iZone], true);
  }


  if(mixingplane){
    if (rank == MASTER_NODE) cout << "Set span-wise sections between zones on Mixing-Plane interface." << endl;
    for (donorZone = 0; donorZone < nZone; donorZone++) {
      for (targetZone = 0; targetZone < nZone; targetZone++) {
        if (targetZone != donorZone){
          interface[donorZone][targetZone]->SetSpanWiseLevels(config[donorZone], config[targetZone]);
        }
      }
    }
  }

  if (rank == MASTER_NODE) cout << "Transfer average geometric quantities to zone 0." << endl;
  for (iZone = 1; iZone < nZone; iZone++) {
    interface[iZone][ZONE_0]->GatherAverageTurboGeoValues(geometry[iZone][INST_0][MESH_0],geometry[ZONE_0][INST_0][MESH_0], iZone);
  }

  /*--- Transfer number of blade to ZONE_0 to correctly compute turbo performance---*/
  for (iZone = 1; iZone < nZone; iZone++) {
    nBlades = config[iZone]->GetnBlades(iZone);
    config[ZONE_0]->SetnBlades(iZone, nBlades);
  }

  if (rank == MASTER_NODE){
    for (iZone = 0; iZone < nZone; iZone++) {
    areaIn  = geometry[iZone][INST_0][MESH_0]->GetSpanAreaIn(iZone, config[iZone]->GetnSpanWiseSections());
    areaOut = geometry[iZone][INST_0][MESH_0]->GetSpanAreaOut(iZone, config[iZone]->GetnSpanWiseSections());
    nBlades = config[iZone]->GetnBlades(iZone);
    cout << "Inlet area for Row "<< iZone + 1<< ": " << areaIn*10000.0 <<" cm^2."  <<endl;
    cout << "Oulet area for Row "<< iZone + 1<< ": " << areaOut*10000.0 <<" cm^2."  <<endl;
    cout << "Recomputed number of blades for Row "<< iZone + 1 << ": " << nBlades<<"."  <<endl;
    }
  }


  if(mixingplane){
    if (rank == MASTER_NODE) cout<<"Preprocessing of the Mixing-Plane Interface." << endl;
    for (donorZone = 0; donorZone < nZone; donorZone++) {
      nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;
      for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
        for (targetZone = 0; targetZone < nZone; targetZone++) {
          if (targetZone != donorZone){
            interface[donorZone][targetZone]->PreprocessAverage(geometry[donorZone][INST_0][MESH_0], geometry[targetZone][INST_0][MESH_0],
                config[donorZone], config[targetZone],
                iMarkerInt);
          }
        }
      }
    }
  }

  if(!restart && !discrete_adjoint){
    if (rank == MASTER_NODE) cout<<"Initialize turbomachinery solution quantities." << endl;
    for(iZone = 0; iZone < nZone; iZone++) {
      solver[iZone][INST_0][MESH_0][FLOW_SOL]->SetFreeStream_TurboSolution(config[iZone]);
    }
  }

  if (rank == MASTER_NODE) cout<<"Initialize inflow and outflow average solution quantities." << endl;
  for(iZone = 0; iZone < nZone; iZone++) {
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->PreprocessAverage(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],INFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->PreprocessAverage(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],OUTFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->TurboAverageProcess(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],INFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->TurboAverageProcess(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],OUTFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->GatherInOutAverageValues(config[iZone], geometry[iZone][INST_0][MESH_0]);
    if (rank == MASTER_NODE){
      flowAngleIn = solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityIn(iZone, config[iZone]->GetnSpanWiseSections())[1];
      flowAngleIn /= solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityIn(iZone, config[iZone]->GetnSpanWiseSections())[0];
      flowAngleIn = atan(flowAngleIn)*180.0/PI_NUMBER;
      cout << "Inlet flow angle for Row "<< iZone + 1<< ": "<< flowAngleIn <<"°."  <<endl;
      flowAngleOut = solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityOut(iZone, config[iZone]->GetnSpanWiseSections())[1];
      flowAngleOut /= solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityOut(iZone, config[iZone]->GetnSpanWiseSections())[0];
      flowAngleOut = atan(flowAngleOut)*180.0/PI_NUMBER;
      cout << "Outlet flow angle for Row "<< iZone + 1<< ": "<< flowAngleOut <<"°."  <<endl;

    }
  }

}

CDriver::~CDriver(void) {}

void CDriver::Print_DirectResidual(RECORDING kind_recording) {

  if (!(rank == MASTER_NODE && kind_recording == RECORDING::SOLUTION_VARIABLES)) return;

  const bool multizone = config_container[ZONE_0]->GetMultizone_Problem();

  /*--- Helper lambda func to return lenghty [iVar][iZone] string.  ---*/
  auto iVar_iZone2string = [&](unsigned short ivar, unsigned short izone) {
    if (multizone)
      return "[" + std::to_string(ivar) + "][" + std::to_string(izone) + "]";
    else
      return "[" + std::to_string(ivar) + "]";
  };

  /*--- Print residuals in the first iteration ---*/

  const unsigned short fieldWidth = 15;
  PrintingToolbox::CTablePrinter RMSTable(&std::cout);
  RMSTable.SetPrecision(config_container[ZONE_0]->GetOutput_Precision());

  /*--- The CTablePrinter requires two sweeps:
    *--- 0. Add the colum names (addVals=0=false) plus CTablePrinter.PrintHeader()
    *--- 1. Add the RMS-residual values (addVals=1=true) plus CTablePrinter.PrintFooter() ---*/
  for (int addVals = 0; addVals < 2; addVals++) {

    for (unsigned short iZone = 0; iZone < nZone; iZone++) {

      auto solvers = solver_container[iZone][INST_0][MESH_0];
      auto configs = config_container[iZone];

      /*--- Note: the FEM-Flow solvers are availalbe for disc. adjoint runs only for SingleZone. ---*/
      if (configs->GetFluidProblem() || configs->GetFEMSolver()) {

        for (unsigned short iVar = 0; iVar < solvers[FLOW_SOL]->GetnVar(); iVar++) {
          if (!addVals)
            RMSTable.AddColumn("rms_Flow" + iVar_iZone2string(iVar, iZone), fieldWidth);
          else
            RMSTable << log10(solvers[FLOW_SOL]->GetRes_RMS(iVar));
        }

        if (configs->GetKind_Turb_Model() != TURB_MODEL::NONE && !configs->GetFrozen_Visc_Disc()) {
          for (unsigned short iVar = 0; iVar < solvers[TURB_SOL]->GetnVar(); iVar++) {
            if (!addVals)
              RMSTable.AddColumn("rms_Turb" + iVar_iZone2string(iVar, iZone), fieldWidth);
            else
              RMSTable << log10(solvers[TURB_SOL]->GetRes_RMS(iVar));
          }
        }

        if (configs->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
          for (unsigned short iVar = 0; iVar < solvers[SPECIES_SOL]->GetnVar(); iVar++) {
            if (!addVals)
              RMSTable.AddColumn("rms_Spec" + iVar_iZone2string(iVar, iZone), fieldWidth);
            else
              RMSTable << log10(solvers[SPECIES_SOL]->GetRes_RMS(iVar));
          }
        }

        if (!multizone && configs->GetWeakly_Coupled_Heat()){
          if (!addVals) RMSTable.AddColumn("rms_Heat" + iVar_iZone2string(0, iZone), fieldWidth);
          else RMSTable << log10(solvers[HEAT_SOL]->GetRes_RMS(0));
        }

        if (configs->AddRadiation()) {
          if (!addVals) RMSTable.AddColumn("rms_Rad" + iVar_iZone2string(0, iZone), fieldWidth);
          else RMSTable << log10(solvers[RAD_SOL]->GetRes_RMS(0));
        }

      }
      else if (configs->GetStructuralProblem()) {

        if (configs->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE){
          if (!addVals) {
            RMSTable.AddColumn("UTOL-A", fieldWidth);
            RMSTable.AddColumn("RTOL-A", fieldWidth);
            RMSTable.AddColumn("ETOL-A", fieldWidth);
          }
          else {
            RMSTable << log10(solvers[FEA_SOL]->GetRes_FEM(0))
                     << log10(solvers[FEA_SOL]->GetRes_FEM(1))
                     << log10(solvers[FEA_SOL]->GetRes_FEM(2));
          }
        }
        else{
          if (!addVals) {
            RMSTable.AddColumn("log10[RMS Ux]", fieldWidth);
            RMSTable.AddColumn("log10[RMS Uy]", fieldWidth);
            if (nDim == 3) RMSTable.AddColumn("log10[RMS Uz]", fieldWidth);
          }
          else {
            RMSTable << log10(solvers[FEA_SOL]->GetRes_FEM(0))
                     << log10(solvers[FEA_SOL]->GetRes_FEM(1));
            if (nDim == 3) RMSTable << log10(solvers[FEA_SOL]->GetRes_FEM(2));
          }
        }

      }
      else if (configs->GetHeatProblem()) {

        if (!addVals) RMSTable.AddColumn("rms_Heat" + iVar_iZone2string(0, iZone), fieldWidth);
        else RMSTable << log10(solvers[HEAT_SOL]->GetRes_RMS(0));
      } else {
        SU2_MPI::Error("Invalid KindSolver for CDiscAdj-MultiZone/SingleZone-Driver.", CURRENT_FUNCTION);
      }
    } // loop iZone

    if (!addVals) RMSTable.PrintHeader();
    else RMSTable.PrintFooter();

  } // for addVals

  cout << "\n-------------------------------------------------------------------------\n" << endl;

}

CFluidDriver::CFluidDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator) : CDriver(confFile, val_nZone, MPICommunicator, false) {
  Max_Iter = config_container[ZONE_0]->GetnInner_Iter();
}

CFluidDriver::~CFluidDriver(void) { }

void CFluidDriver::StartSolver(){

#ifdef VTUNEPROF
  __itt_resume();
#endif

  /*--- Main external loop of the solver. Within this loop, each iteration ---*/

  if (rank == MASTER_NODE){
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
  }

  unsigned long Iter = 0;
  while ( Iter < Max_Iter ) {

    /*--- Perform some external iteration preprocessing. ---*/

    Preprocess(Iter);

    /*--- Perform a dynamic mesh update if required. ---*/
    /*--- For the Disc.Adj. of a case with (rigidly) moving grid, the appropriate
          mesh cordinates are read from the restart files. ---*/
    if (!fem_solver &&
        !(config_container[ZONE_0]->GetGrid_Movement() && config_container[ZONE_0]->GetDiscrete_Adjoint())) {
      DynamicMeshUpdate(Iter);
    }

    /*--- Run a single iteration of the problem (fluid, elasticity, heat, ...). ---*/

    Run();

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Terminate the simulation if only the Jacobian must be computed. ---*/
    if (config_container[ZONE_0]->GetJacobian_Spatial_Discretization_Only()) break;

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(Iter);

    /*--- Output the solution in files. ---*/

    Output(Iter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    Iter++;

  }
#ifdef VTUNEPROF
  __itt_pause();
#endif
}


void CFluidDriver::Preprocess(unsigned long Iter) {

  /*--- Set the value of the external iteration and physical time. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetInnerIter(Iter);
    if (config_container[iZone]->GetTime_Marching() != TIME_MARCHING::STEADY)
      config_container[iZone]->SetPhysicalTime(static_cast<su2double>(Iter)*config_container[iZone]->GetDelta_UnstTimeND());
    else
      config_container[iZone]->SetPhysicalTime(0.0);
  }

  /*--- Set the initial condition for EULER/N-S/RANS and for a non FSI simulation ---*/

  if(!fsi) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone]->GetFluidProblem()) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++) {
          solver_container[iZone][iInst][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone][INST_0], solver_container[iZone][iInst], config_container[iZone], Iter);
        }
      }
    }
  }
}

void CFluidDriver::Run() {

  unsigned short iZone, jZone, checkConvergence;
  unsigned long IntIter, nIntIter;
  bool unsteady;

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
             (config_container[MESH_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

  /*--- Updating zone interface communication patterns,
   needed only for unsteady simulation since for steady problems
   this is done once in the interpolator_container constructor
   at the beginning of the computation ---*/

  if ( unsteady ) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && interpolator_container[iZone][jZone] != nullptr)
        interpolator_container[iZone][jZone]->SetTransferCoeff(config_container);
    }
  }

  /*--- Begin Unsteady pseudo-time stepping internal loop, if not unsteady it does only one step --*/

  if (unsteady)
    nIntIter = config_container[MESH_0]->GetnInner_Iter();
  else
    nIntIter = 1;

  for (IntIter = 0; IntIter < nIntIter; IntIter++) {

    /*--- At each pseudo time-step updates transfer data ---*/
    for (iZone = 0; iZone < nZone; iZone++)
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && interface_container[iZone][jZone] != nullptr)
          Transfer_Data(iZone, jZone);

    /*--- For each zone runs one single iteration ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetInnerIter(IntIter);
      iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container, solver_container, numerics_container,
                                                  config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
    }

    /*--- Check convergence in each zone --*/

    checkConvergence = 0;
    for (iZone = 0; iZone < nZone; iZone++)
    checkConvergence += (int) integration_container[iZone][INST_0][FLOW_SOL]->GetConvergence();

    /*--- If convergence was reached in every zone --*/

  if (checkConvergence == nZone) break;
  }

}

void CFluidDriver::Transfer_Data(unsigned short donorZone, unsigned short targetZone) {

  interface_container[donorZone][targetZone]->BroadcastData(*interpolator_container[donorZone][targetZone].get(),
    solver_container[donorZone][INST_0][MESH_0][FLOW_SOL], solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
    geometry_container[donorZone][INST_0][MESH_0], geometry_container[targetZone][INST_0][MESH_0],
    config_container[donorZone], config_container[targetZone]);

  if (config_container[targetZone]->GetKind_Solver() == MAIN_SOLVER::RANS) {
    interface_container[donorZone][targetZone]->BroadcastData(*interpolator_container[donorZone][targetZone].get(),
      solver_container[donorZone][INST_0][MESH_0][TURB_SOL], solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
      geometry_container[donorZone][INST_0][MESH_0], geometry_container[targetZone][INST_0][MESH_0],
      config_container[donorZone], config_container[targetZone]);
  }
}

void CFluidDriver::Update() {

  for(iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone][INST_0]->Update(output_container[iZone], integration_container, geometry_container,
         solver_container, numerics_container, config_container,
         surface_movement, grid_movement, FFDBox, iZone, INST_0);
}

void CFluidDriver::DynamicMeshUpdate(unsigned long TimeIter) {

  bool harmonic_balance;

  for (iZone = 0; iZone < nZone; iZone++) {
   harmonic_balance = (config_container[iZone]->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance)) {
      iteration_container[iZone][INST_0]->SetGrid_Movement(geometry_container[iZone][INST_0], surface_movement[iZone], grid_movement[iZone][INST_0], solver_container[iZone][INST_0], config_container[iZone], 0, TimeIter );
    }
  }

}
bool CFluidDriver::Monitor(unsigned long ExtIter) {

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  IterCount++;
  UsedTime = (StopTime - StartTime) + UsedTimeCompute;

  /*--- Check if there is any change in the runtime parameters ---*/

  CConfig *runtime = nullptr;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetTimeIter(ExtIter);
  delete runtime;

  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
    case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::NEMO_EULER: case MAIN_SOLVER::NEMO_NAVIER_STOKES:
      StopCalc = integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence(); break;
    case MAIN_SOLVER::HEAT_EQUATION:
      StopCalc = integration_container[ZONE_0][INST_0][HEAT_SOL]->GetConvergence(); break;
    case MAIN_SOLVER::FEM_ELASTICITY:
      StopCalc = integration_container[ZONE_0][INST_0][FEA_SOL]->GetConvergence(); break;
    case MAIN_SOLVER::ADJ_EULER: case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER: case MAIN_SOLVER::DISC_ADJ_FEM_NS: case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
      StopCalc = integration_container[ZONE_0][INST_0][ADJFLOW_SOL]->GetConvergence(); break;
    default:
      break;
  }

  /*--- Set StopCalc to true if max. number of iterations has been reached ---*/

  StopCalc = StopCalc || (ExtIter == Max_Iter - 1);

  return StopCalc;

}


void CFluidDriver::Output(unsigned long InnerIter) {

  for (iZone = 0; iZone < nZone; iZone++) {
    const auto inst = config_container[iZone]->GetiInst();

    for (iInst = 0; iInst < nInst[iZone]; ++iInst) {
      config_container[iZone]->SetiInst(iInst);
      output_container[iZone]->SetResult_Files(geometry_container[iZone][iInst][MESH_0],
                                               config_container[iZone],
                                               solver_container[iZone][iInst][MESH_0],
                                               InnerIter, StopCalc);
    }
    config_container[iZone]->SetiInst(inst);
  }

}


CTurbomachineryDriver::CTurbomachineryDriver(char* confFile, unsigned short val_nZone,
                                             SU2_Comm MPICommunicator):
                                             CFluidDriver(confFile, val_nZone, MPICommunicator) {

  output_legacy = COutputFactory::CreateLegacyOutput(config_container[ZONE_0]);

  /*--- LEGACY OUTPUT (going to be removed soon) --- */

  /*--- Open the convergence history file ---*/
  ConvHist_file = nullptr;
  ConvHist_file = new ofstream*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    ConvHist_file[iZone] = nullptr;
    if (rank == MASTER_NODE){
      ConvHist_file[iZone] = new ofstream[nInst[iZone]];
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        output_legacy->SetConvHistory_Header(&ConvHist_file[iZone][iInst], config_container[iZone], iZone, iInst);
      }
    }
  }

  if (nZone > 1){
    Max_Iter = config_container[ZONE_0]->GetnOuter_Iter();
  }
}

CTurbomachineryDriver::~CTurbomachineryDriver(void) {
  if (rank == MASTER_NODE){
    /*--- Close the convergence history file. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < 1; iInst++) {
        ConvHist_file[iZone][iInst].close();
      }
      delete [] ConvHist_file[iZone];
    }
    delete [] ConvHist_file;
  }
}

void CTurbomachineryDriver::Run() {

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  /* --- Update the mixing-plane interface ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if(mixingplane)SetMixingPlane(iZone);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Postprocess(output_container[iZone], integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  if (rank == MASTER_NODE){
    SetTurboPerformance(ZONE_0);
  }


}

void CTurbomachineryDriver::SetMixingPlane(unsigned short donorZone){

  unsigned short targetZone, nMarkerInt, iMarkerInt ;
  nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;

  /* --- transfer the average value from the donorZone to the targetZone*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
    for (targetZone = 0; targetZone < nZone; targetZone++) {
      if (targetZone != donorZone){
        interface_container[donorZone][targetZone]->AllgatherAverage(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone], iMarkerInt );
      }
    }
  }
}

void CTurbomachineryDriver::SetTurboPerformance(unsigned short targetZone){

  unsigned short donorZone;
  //IMPORTANT this approach of multi-zone performances rely upon the fact that turbomachinery markers follow the natural (stator-rotor) development of the real machine.
  /* --- transfer the local turboperfomance quantities (for each blade)  from all the donorZones to the targetZone (ZONE_0) ---*/
  for (donorZone = 1; donorZone < nZone; donorZone++) {
    interface_container[donorZone][targetZone]->GatherAverageValues(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL], donorZone);
  }

  /* --- compute turboperformance for each stage and the global machine ---*/

 output_legacy->ComputeTurboPerformance(solver_container[targetZone][INST_0][MESH_0][FLOW_SOL], geometry_container[targetZone][INST_0][MESH_0], config_container[targetZone]);

}


bool CTurbomachineryDriver::Monitor(unsigned long ExtIter) {

  su2double rot_z_ini, rot_z_final ,rot_z;
  su2double outPres_ini, outPres_final, outPres;
  unsigned long rampFreq, finalRamp_Iter;
  unsigned short iMarker, KindBC, KindBCOption;
  string Marker_Tag;

  bool print;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  IterCount++;
  UsedTime = (StopTime - StartTime);


  /*--- Check if there is any change in the runtime parameters ---*/
  CConfig *runtime = nullptr;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetInnerIter(ExtIter);
  delete runtime;

  /*--- Update the convergence history file (serial and parallel computations). ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      output_legacy->SetConvHistory_Body(&ConvHist_file[iZone][iInst], geometry_container, solver_container,
          config_container, integration_container, false, UsedTime, iZone, iInst);
  }

  /*--- ROTATING FRAME Ramp: Compute the updated rotational velocity. ---*/
  if (config_container[ZONE_0]->GetGrid_Movement() && config_container[ZONE_0]->GetRampRotatingFrame()) {
    rampFreq       = SU2_TYPE::Int(config_container[ZONE_0]->GetRampRotatingFrame_Coeff(1));
    finalRamp_Iter = SU2_TYPE::Int(config_container[ZONE_0]->GetRampRotatingFrame_Coeff(2));
    rot_z_ini = config_container[ZONE_0]->GetRampRotatingFrame_Coeff(0);
    print = false;
    if(ExtIter % rampFreq == 0 &&  ExtIter <= finalRamp_Iter){

      for (iZone = 0; iZone < nZone; iZone++) {
        rot_z_final = config_container[iZone]->GetFinalRotation_Rate_Z();
        if(abs(rot_z_final) > 0.0){
          rot_z = rot_z_ini + ExtIter*( rot_z_final - rot_z_ini)/finalRamp_Iter;
          config_container[iZone]->SetRotation_Rate(2, rot_z);
          if(rank == MASTER_NODE && print && ExtIter > 0) {
            cout << endl << " Updated rotating frame grid velocities";
            cout << " for zone " << iZone << "." << endl;
          }
          geometry_container[iZone][INST_0][MESH_0]->SetRotationalVelocity(config_container[iZone], print);
          geometry_container[iZone][INST_0][MESH_0]->SetShroudVelocity(config_container[iZone]);
        }
      }

      for (iZone = 0; iZone < nZone; iZone++) {
        geometry_container[iZone][INST_0][MESH_0]->SetAvgTurboValue(config_container[iZone], iZone, INFLOW, false);
        geometry_container[iZone][INST_0][MESH_0]->SetAvgTurboValue(config_container[iZone],iZone, OUTFLOW, false);
        geometry_container[iZone][INST_0][MESH_0]->GatherInOutAverageValues(config_container[iZone], false);

      }

      for (iZone = 1; iZone < nZone; iZone++) {
        interface_container[iZone][ZONE_0]->GatherAverageTurboGeoValues(geometry_container[iZone][INST_0][MESH_0],geometry_container[ZONE_0][INST_0][MESH_0], iZone);
      }

    }
  }


  /*--- Outlet Pressure Ramp: Compute the updated rotational velocity. ---*/
  if (config_container[ZONE_0]->GetRampOutletPressure()) {
    rampFreq       = SU2_TYPE::Int(config_container[ZONE_0]->GetRampOutletPressure_Coeff(1));
    finalRamp_Iter = SU2_TYPE::Int(config_container[ZONE_0]->GetRampOutletPressure_Coeff(2));
    outPres_ini    = config_container[ZONE_0]->GetRampOutletPressure_Coeff(0);
    outPres_final  = config_container[ZONE_0]->GetFinalOutletPressure();

    if(ExtIter % rampFreq == 0 &&  ExtIter <= finalRamp_Iter){
      outPres = outPres_ini + ExtIter*(outPres_final - outPres_ini)/finalRamp_Iter;
      if(rank == MASTER_NODE) config_container[ZONE_0]->SetMonitotOutletPressure(outPres);

      for (iZone = 0; iZone < nZone; iZone++) {
        for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
          KindBC = config_container[iZone]->GetMarker_All_KindBC(iMarker);
          switch (KindBC) {
          case RIEMANN_BOUNDARY:
            Marker_Tag         = config_container[iZone]->GetMarker_All_TagBound(iMarker);
            KindBCOption       = config_container[iZone]->GetKind_Data_Riemann(Marker_Tag);
            if(KindBCOption == STATIC_PRESSURE || KindBCOption == RADIAL_EQUILIBRIUM ){
              SU2_MPI::Error("Outlet pressure ramp only implemented for NRBC", CURRENT_FUNCTION);
            }
            break;
          case GILES_BOUNDARY:
            Marker_Tag         = config_container[iZone]->GetMarker_All_TagBound(iMarker);
            KindBCOption       = config_container[iZone]->GetKind_Data_Giles(Marker_Tag);
            if(KindBCOption == STATIC_PRESSURE || KindBCOption == STATIC_PRESSURE_1D || KindBCOption == RADIAL_EQUILIBRIUM ){
              config_container[iZone]->SetGiles_Var1(outPres, Marker_Tag);
            }
            break;
          }
        }
      }
    }
  }


  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
  case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
  case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
  case MAIN_SOLVER::NEMO_EULER: case MAIN_SOLVER::NEMO_NAVIER_STOKES:
    StopCalc = integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence();
    break;
  case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
  case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
  case MAIN_SOLVER::DISC_ADJ_FEM_EULER: case MAIN_SOLVER::DISC_ADJ_FEM_NS: case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
    StopCalc = integration_container[ZONE_0][INST_0][ADJFLOW_SOL]->GetConvergence();
    break;
  default:
    break;

  }

  /*--- Set StopCalc to true if max. number of iterations has been reached ---*/

  StopCalc = StopCalc || (ExtIter == Max_Iter - 1);

  return StopCalc;

}

CHBDriver::CHBDriver(char* confFile,
    unsigned short val_nZone,
    SU2_Comm MPICommunicator) : CFluidDriver(confFile,
        val_nZone,
        MPICommunicator) {
  unsigned short kInst;

  nInstHB = nInst[ZONE_0];

  D = nullptr;
  /*--- allocate dynamic memory for the Harmonic Balance operator ---*/
  D = new su2double*[nInstHB]; for (kInst = 0; kInst < nInstHB; kInst++) D[kInst] = new su2double[nInstHB];

  output_legacy = COutputFactory::CreateLegacyOutput(config_container[ZONE_0]);

  /*--- Open the convergence history file ---*/
  ConvHist_file = nullptr;
  ConvHist_file = new ofstream*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    ConvHist_file[iZone] = nullptr;
    if (rank == MASTER_NODE){
      ConvHist_file[iZone] = new ofstream[nInst[iZone]];
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        output_legacy->SetConvHistory_Header(&ConvHist_file[iZone][iInst], config_container[iZone], iZone, iInst);
      }
    }
  }


}

CHBDriver::~CHBDriver(void) {

  unsigned short kInst;

  /*--- delete dynamic memory for the Harmonic Balance operator ---*/
  for (kInst = 0; kInst < nInstHB; kInst++) delete [] D[kInst];
  delete [] D;

  if (rank == MASTER_NODE){
  /*--- Close the convergence history file. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInstHB; iInst++) {
      ConvHist_file[iZone][iInst].close();
    }
    delete [] ConvHist_file[iZone];
  }
  delete [] ConvHist_file;
  }
}


void CHBDriver::Run() {

  /*--- Run a single iteration of a Harmonic Balance problem. Preprocess all
   all zones before beginning the iteration. ---*/

  for (iInst = 0; iInst < nInstHB; iInst++)
    iteration_container[ZONE_0][iInst]->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  for (iInst = 0; iInst < nInstHB; iInst++)
    iteration_container[ZONE_0][iInst]->Iterate(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  /*--- Update the convergence history file (serial and parallel computations). ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      output_legacy->SetConvHistory_Body(&ConvHist_file[iZone][iInst], geometry_container, solver_container,
          config_container, integration_container, false, UsedTime, iZone, iInst);
  }

}

void CHBDriver::Update() {

  for (iInst = 0; iInst < nInstHB; iInst++) {
    /*--- Compute the harmonic balance terms across all zones ---*/
    SetHarmonicBalance(iInst);

  }

  /*--- Precondition the harmonic balance source terms ---*/
  if (config_container[ZONE_0]->GetHB_Precondition() == YES) {
    StabilizeHarmonicBalance();

  }

  for (iInst = 0; iInst < nInstHB; iInst++) {

    /*--- Update the harmonic balance terms across all zones ---*/
    iteration_container[ZONE_0][iInst]->Update(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  }

}

void CHBDriver::ResetConvergence() {

  for(iInst = 0; iInst < nZone; iInst++) {
    switch (config_container[ZONE_0]->GetKind_Solver()) {

    case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
      integration_container[ZONE_0][iInst][FLOW_SOL]->SetConvergence(false);
      if (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::RANS) integration_container[ZONE_0][iInst][TURB_SOL]->SetConvergence(false);
      if(config_container[ZONE_0]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) integration_container[ZONE_0][iInst][TRANS_SOL]->SetConvergence(false);
      break;

    case MAIN_SOLVER::FEM_ELASTICITY:
      integration_container[ZONE_0][iInst][FEA_SOL]->SetConvergence(false);
      break;

    case MAIN_SOLVER::ADJ_EULER: case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
      integration_container[ZONE_0][iInst][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) || (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) )
        integration_container[ZONE_0][iInst][ADJTURB_SOL]->SetConvergence(false);
      break;

    default:
      SU2_MPI::Error("Harmonic Balance has not been set up for this solver.", CURRENT_FUNCTION);
    }
  }

}

void CHBDriver::SetHarmonicBalance(unsigned short iInst) {

  unsigned short iVar, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
  if (adjoint) {
    implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }

  unsigned long InnerIter = config_container[ZONE_0]->GetInnerIter();

  /*--- Retrieve values from the config file ---*/
  su2double *U = new su2double[nVar];
  su2double *U_old = new su2double[nVar];
  su2double *Psi = new su2double[nVar];
  su2double *Psi_old = new su2double[nVar];
  su2double *Source = new su2double[nVar];
  su2double deltaU, deltaPsi;

  /*--- Compute period of oscillation ---*/
  su2double period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.  */
  period /= config_container[ZONE_0]->GetTime_Ref();

  if (InnerIter == 0)
    ComputeHB_Operator();

  /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][iInst][iMGlevel]->GetnPoint(); iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {
        Source[iVar] = 0.0;
      }

      /*--- Step across the columns ---*/
      for (jInst = 0; jInst < nInstHB; jInst++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar; iVar++) {

          if (!adjoint) {
            U[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
            Source[iVar] += U[iVar]*D[iInst][jInst];

            if (implicit) {
              U_old[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
              deltaU = U[iVar] - U_old[iVar];
              Source[iVar] += deltaU*D[iInst][jInst];
            }

          }

          else {
            Psi[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
            Source[iVar] += Psi[iVar]*D[jInst][iInst];

            if (implicit) {
              Psi_old[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
              deltaPsi = Psi[iVar] - Psi_old[iVar];
              Source[iVar] += deltaPsi*D[jInst][iInst];
            }
          }
        }

        /*--- Store sources for current row ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (!adjoint) {
            solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iVar]);
          }
          else {
            solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iVar]);
          }
        }

      }
    }
  }

  /*--- Source term for a turbulence model ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::RANS) {

    /*--- Extra variables needed if we have a turbulence model. ---*/
    unsigned short nVar_Turb = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar();
    su2double *U_Turb = new su2double[nVar_Turb];
    su2double *Source_Turb = new su2double[nVar_Turb];

    /*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
      for (jInst = 0; jInst < nInstHB; jInst++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar_Turb; iVar++) {
          U_Turb[iVar] = solver_container[ZONE_0][jInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint, iVar);
          Source_Turb[iVar] += U_Turb[iVar]*D[iInst][jInst];
        }
      }

      /*--- Store sources for current iZone ---*/
      for (iVar = 0; iVar < nVar_Turb; iVar++)
        solver_container[ZONE_0][iInst][MESH_0][TURB_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source_Turb[iVar]);
    }

    delete [] U_Turb;
    delete [] Source_Turb;
  }

  delete [] Source;
  delete [] U;
  delete [] U_old;
  delete [] Psi;
  delete [] Psi_old;

}

void CHBDriver::StabilizeHarmonicBalance() {

  unsigned short i, j, k, iVar, iInst, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());

  /*--- Retrieve values from the config file ---*/
  su2double *Source     = new su2double[nInstHB];
  su2double *Source_old = new su2double[nInstHB];
  su2double Delta;

  su2double **Pinv     = new su2double*[nInstHB];
  su2double **P        = new su2double*[nInstHB];
  for (iInst = 0; iInst < nInstHB; iInst++) {
    Pinv[iInst]       = new su2double[nInstHB];
    P[iInst]          = new su2double[nInstHB];
  }

  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][iMGlevel]->GetnPoint(); iPoint++) {

      /*--- Get time step for current node ---*/
      Delta = solver_container[ZONE_0][INST_0][iMGlevel][FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint);

      /*--- Setup stabilization matrix for this node ---*/
      for (iInst = 0; iInst < nInstHB; iInst++) {
        for (jInst = 0; jInst < nInstHB; jInst++) {
          if (jInst == iInst ) {
            Pinv[iInst][jInst] = 1.0 + Delta*D[iInst][jInst];
          }
          else {
            Pinv[iInst][jInst] = Delta*D[iInst][jInst];
          }
        }
      }

      /*--- Invert stabilization matrix Pinv with Gauss elimination---*/

      /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
      su2double **temp = new su2double*[nInstHB];
      for (i = 0; i < nInstHB; i++) {
        temp[i] = new su2double[2 * nInstHB];
      }

      /*---  Copy the desired matrix into the temporary matrix ---*/
      for (i = 0; i < nInstHB; i++) {
        for (j = 0; j < nInstHB; j++) {
          temp[i][j] = Pinv[i][j];
          temp[i][nInstHB + j] = 0;
        }
        temp[i][nInstHB + i] = 1;
      }

      su2double max_val;
      unsigned short max_idx;

      /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
      for (k = 0; k < nInstHB - 1; k++) {
        max_idx = k;
        max_val = abs(temp[k][k]);
        /*---  Find the largest value (pivot) in the column  ---*/
        for (j = k; j < nInstHB; j++) {
          if (abs(temp[j][k]) > max_val) {
            max_idx = j;
            max_val = abs(temp[j][k]);
          }
        }

        /*---  Move the row with the highest value up  ---*/
        for (j = 0; j < (nInstHB * 2); j++) {
          su2double d = temp[k][j];
          temp[k][j] = temp[max_idx][j];
          temp[max_idx][j] = d;
        }
        /*---  Subtract the moved row from all other rows ---*/
        for (i = k + 1; i < nInstHB; i++) {
          su2double c = temp[i][k] / temp[k][k];
          for (j = 0; j < (nInstHB * 2); j++) {
            temp[i][j] = temp[i][j] - temp[k][j] * c;
          }
        }
      }

      /*---  Back-substitution  ---*/
      for (k = nInstHB - 1; k > 0; k--) {
        if (temp[k][k] != su2double(0.0)) {
          for (int i = k - 1; i > -1; i--) {
            su2double c = temp[i][k] / temp[k][k];
            for (j = 0; j < (nInstHB * 2); j++) {
              temp[i][j] = temp[i][j] - temp[k][j] * c;
            }
          }
        }
      }

      /*---  Normalize the inverse  ---*/
      for (i = 0; i < nInstHB; i++) {
        su2double c = temp[i][i];
        for (j = 0; j < nInstHB; j++) {
          temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
        }
      }

      /*---  Copy the inverse back to the main program flow ---*/
      for (i = 0; i < nInstHB; i++) {
        for (j = 0; j < nInstHB; j++) {
          P[i][j] = temp[i][j + nInstHB];
        }
      }

      /*---  Delete dynamic template  ---*/
      for (iInst = 0; iInst < nInstHB; iInst++) {
        delete[] temp[iInst];
      }
      delete[] temp;

      /*--- Loop through variables to precondition ---*/
      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Get current source terms (not yet preconditioned) and zero source array to prepare preconditioning ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          Source_old[iInst] = solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->GetHarmonicBalance_Source(iPoint, iVar);
          Source[iInst] = 0;
        }

        /*--- Step through columns ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          for (jInst = 0; jInst < nInstHB; jInst++) {
            Source[iInst] += P[iInst][jInst]*Source_old[jInst];
          }

          /*--- Store updated source terms for current node ---*/
          if (!adjoint) {
            solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
          }
          else {
            solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
          }
        }

      }
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iInst = 0; iInst < nInstHB; iInst++){
    delete [] P[iInst];
    delete [] Pinv[iInst];
  }
  delete [] P;
  delete [] Pinv;
  delete [] Source;
  delete [] Source_old;

}

void CHBDriver::ComputeHB_Operator() {

  const   complex<su2double> J(0.0,1.0);
  unsigned short i, j, k, iInst;

  su2double *Omega_HB       = new su2double[nInstHB];
  complex<su2double> **E    = new complex<su2double>*[nInstHB];
  complex<su2double> **Einv = new complex<su2double>*[nInstHB];
  complex<su2double> **DD   = new complex<su2double>*[nInstHB];
  for (iInst = 0; iInst < nInstHB; iInst++) {
    E[iInst]    = new complex<su2double>[nInstHB];
    Einv[iInst] = new complex<su2double>[nInstHB];
    DD[iInst]   = new complex<su2double>[nInstHB];
  }

  /*--- Get simualation period from config file ---*/
  su2double Period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.      */
  Period /= config_container[ZONE_0]->GetTime_Ref();

  /*--- Build the array containing the selected frequencies to solve ---*/
  for (iInst = 0; iInst < nInstHB; iInst++) {
    Omega_HB[iInst]  = config_container[ZONE_0]->GetOmega_HB()[iInst];
    Omega_HB[iInst] /= config_container[ZONE_0]->GetOmega_Ref(); //TODO: check
  }

  /*--- Build the diagonal matrix of the frequencies DD ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      if (k == i ) {
        DD[i][k] = J*Omega_HB[k];
      }
    }
  }


  /*--- Build the harmonic balance inverse matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      Einv[i][k] = complex<su2double>(cos(Omega_HB[k]*(i*Period/nInstHB))) + J*complex<su2double>(sin(Omega_HB[k]*(i*Period/nInstHB)));
    }
  }

  /*---  Invert inverse harmonic balance Einv with Gauss elimination ---*/

  /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
  complex<su2double> **temp = new complex<su2double>*[nInstHB];
  for (i = 0; i < nInstHB; i++) {
    temp[i] = new complex<su2double>[2 * nInstHB];
  }

  /*---  Copy the desired matrix into the temporary matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (j = 0; j < nInstHB; j++) {
      temp[i][j] = Einv[i][j];
      temp[i][nInstHB + j] = 0;
    }
    temp[i][nInstHB + i] = 1;
  }

  su2double max_val;
  unsigned short max_idx;

  /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
  for (k = 0; k < nInstHB - 1; k++) {
    max_idx = k;
    max_val = abs(temp[k][k]);
    /*---  Find the largest value (pivot) in the column  ---*/
    for (j = k; j < nInstHB; j++) {
      if (abs(temp[j][k]) > max_val) {
        max_idx = j;
        max_val = abs(temp[j][k]);
      }
    }
    /*---  Move the row with the highest value up  ---*/
    for (j = 0; j < (nInstHB * 2); j++) {
      complex<su2double> d = temp[k][j];
      temp[k][j] = temp[max_idx][j];
      temp[max_idx][j] = d;
    }
    /*---  Subtract the moved row from all other rows ---*/
    for (i = k + 1; i < nInstHB; i++) {
      complex<su2double> c = temp[i][k] / temp[k][k];
      for (j = 0; j < (nInstHB * 2); j++) {
        temp[i][j] = temp[i][j] - temp[k][j] * c;
      }
    }
  }
  /*---  Back-substitution  ---*/
  for (k = nInstHB - 1; k > 0; k--) {
    if (temp[k][k] != complex<su2double>(0.0)) {
      for (int i = k - 1; i > -1; i--) {
        complex<su2double> c = temp[i][k] / temp[k][k];
        for (j = 0; j < (nInstHB * 2); j++) {
          temp[i][j] = temp[i][j] - temp[k][j] * c;
        }
      }
    }
  }
  /*---  Normalize the inverse  ---*/
  for (i = 0; i < nInstHB; i++) {
    complex<su2double> c = temp[i][i];
    for (j = 0; j < nInstHB; j++) {
      temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
    }
  }
  /*---  Copy the inverse back to the main program flow ---*/
  for (i = 0; i < nInstHB; i++) {
    for (j = 0; j < nInstHB; j++) {
      E[i][j] = temp[i][j + nInstHB];
    }
  }
  /*---  Delete dynamic template  ---*/
  for (i = 0; i < nInstHB; i++) {
    delete[] temp[i];
  }
  delete[] temp;


  /*---  Temporary matrix for performing product  ---*/
  complex<su2double> **Temp    = new complex<su2double>*[nInstHB];

  /*---  Temporary complex HB operator  ---*/
  complex<su2double> **Dcpx    = new complex<su2double>*[nInstHB];

  for (iInst = 0; iInst < nInstHB; iInst++){
    Temp[iInst]    = new complex<su2double>[nInstHB];
    Dcpx[iInst]   = new complex<su2double>[nInstHB];
  }


  /*---  Calculation of the HB operator matrix ---*/
  for (int row = 0; row < nInstHB; row++) {
    for (int col = 0; col < nInstHB; col++) {
      for (int inner = 0; inner < nInstHB; inner++) {
        Temp[row][col] += Einv[row][inner] * DD[inner][col];
      }
    }
  }

  unsigned short row, col, inner;

  for (row = 0; row < nInstHB; row++) {
    for (col = 0; col < nInstHB; col++) {
      for (inner = 0; inner < nInstHB; inner++) {
        Dcpx[row][col] += Temp[row][inner] * E[inner][col];
      }
    }
  }

  /*---  Take just the real part of the HB operator matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      D[i][k] = real(Dcpx[i][k]);
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iInst = 0; iInst < nInstHB; iInst++){
    delete [] E[iInst];
    delete [] Einv[iInst];
    delete [] DD[iInst];
    delete [] Temp[iInst];
    delete [] Dcpx[iInst];
  }
  delete [] E;
  delete [] Einv;
  delete [] DD;
  delete [] Temp;
  delete [] Dcpx;
  delete [] Omega_HB;

}
