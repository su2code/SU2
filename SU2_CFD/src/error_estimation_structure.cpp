/*!
 * \file error_estimation_structure.cpp
 * \brief Functions for error estimation.
 * \author B. Munguía
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

#include "../include/error_estimation_structure.hpp"

CErrorEstimationDriver::CErrorEstimationDriver(char* confFile,
                                               unsigned short val_nZone,
                                               unsigned short val_nDim, 
                                               bool val_periodic,
                                               SU2_Comm MPICommunicator):config_file_name(confFile), nZone(val_nZone), nDim(val_nDim), fsi(false), fem_solver(false) {

  char zone_file_name[MAX_STRING_SIZE];

  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  geometry_container            = NULL;
  solver_container              = NULL;
  config_container              = NULL;
  output                        = NULL;
  nInst                         = NULL;

  /*--- Definition of the containers for all possible zones. ---*/

  geometry_container            = new CGeometry***[nZone];
  solver_container              = new CSolver****[nZone];
  config_container              = new CConfig*[nZone];
  nInst                         = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]                   = NULL;
    geometry_container[iZone]                 = NULL;
    config_container[iZone]                   = NULL;
    nInst[iZone]                              = 1;
  }

  /*--- Initialize the configuration of the driver ---*/

  driver_config = new CConfig(config_file_name, SU2_MET, ZONE_0, nZone, nDim, VERB_NONE);

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (driver_config->GetKind_Solver() == MULTIZONE){
      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config_container[iZone] = new CConfig(zone_file_name, SU2_MET, iZone, nZone, nDim, VERB_HIGH);
    }
    else{
      config_container[iZone] = new CConfig(config_file_name, SU2_MET, iZone, nZone, nDim, VERB_HIGH);
    }

    /*--- Set the MPI communicator ---*/

    config_container[iZone]->SetMPICommunicator(MPICommunicator);

    /* --- For the output config, enable restart reading (to export sensor fields) and change grid file to target mesh ---*/

    config_container[iZone]->SetRestart(true);

    config_container[iZone]->SetMGLevels(0);

  }

  /*--- Set the multizone part of the problem. ---*/
  if (driver_config->GetKind_Solver() == MULTIZONE){
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }

  /*--- Preprocessing of the config and mesh files. In this routine, the config file is read
   and it is determined whether a problem is single physics or multiphysics. . ---*/

  Input_Preprocessing(config_container, geometry_container, val_periodic);

  /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
   based data structure is constructed, i.e. node and cell neighbors are
   identified and linked, and face areas are computed ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;

  /*--- Determine whether or not the FEM solver is used, which decides the
   type of geometry classes that are instantiated. Only adapted for single-zone problems ---*/
  fem_solver = ((config_container[ZONE_0]->GetKind_Solver() == FEM_EULER)          ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_NAVIER_STOKES)  ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_RANS)           ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_LES)            ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

  if( fem_solver ) {
    switch( config_container[ZONE_0]->GetKind_FEM_Flow() ) {
      case DG: {
        Geometrical_Preprocessing_DGFEM(config_container, geometry_container);
        break;
      }
    }
  }
  else {
    Geometrical_Preprocessing(config_container, geometry_container);
  }

  for (iZone = 0; iZone < nZone; iZone++) {

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Computation of positive surface area in the z-plane which is used for
        the calculation of force coefficient (non-dimensionalization). ---*/

      geometry_container[iZone][iInst][MESH_0]->SetPositive_ZArea(config_container[iZone]);

    }

  }

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Definition of the solver class: solver_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS].
    The solver classes are specific to a particular set of governing equations,
    and they contain the subroutines with instructions for computing each spatial
    term of the PDE, i.e. loops over the edges to compute convective and viscous
    fluxes, loops over the nodes to compute source terms, and routines for
    imposing various boundary condition type for the PDE. ---*/

    if (rank == MASTER_NODE)
      cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;

    solver_container[iZone] = new CSolver*** [nInst[iZone]];


    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      solver_container[iZone][iInst] = NULL;
      solver_container[iZone][iInst] = new CSolver** [MESH_0+1];
      solver_container[iZone][iInst][MESH_0] = NULL;
      solver_container[iZone][iInst][MESH_0] = new CSolver* [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        solver_container[iZone][iInst][MESH_0][iSol] = NULL;
      
      Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
                           config_container[iZone], iInst);

    } // End of loop over iInst

  }

  /*--- Definition of the output class (one for all zones). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  output = new COutput(config_container[ZONE_0]);

}

void CErrorEstimationDriver::Input_Preprocessing(CConfig **config_container, CGeometry ****geometry_container, bool val_periodic) {

  bool fem_solver = false; 

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Determine whether or not the FEM solver is used, which decides the
     type of geometry classes that are instantiated. ---*/
    fem_solver = ((config_container[iZone]->GetKind_Solver() == FEM_EULER)         ||
                  (config_container[iZone]->GetKind_Solver() == FEM_NAVIER_STOKES) ||
                  (config_container[iZone]->GetKind_Solver() == FEM_RANS)          ||
                  (config_container[iZone]->GetKind_Solver() == FEM_LES)           ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone] = new CGeometry** [nInst[iZone]];

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      config_container[iZone]->SetiInst(iInst);

      /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/

      CGeometry *geometry_aux = NULL;

      /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

      geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

      /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

      if ( fem_solver ) geometry_aux->SetColorFEMGrid_Parallel(config_container[iZone]);
      else              geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

      /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

      geometry_container[iZone][iInst]            = NULL;
      geometry_container[iZone][iInst]            = new CGeometry *[MESH_0+1];
      geometry_container[iZone][iInst][MESH_0]    = NULL;


      if( fem_solver ) {
        switch( config_container[iZone]->GetKind_FEM_Flow() ) {
          case DG: {
            geometry_container[iZone][iInst][MESH_0] = new CMeshFEM_DG(geometry_aux, config_container[iZone]);
            break;
          }

          default: {
            SU2_MPI::Error("Unknown FEM flow solver.", CURRENT_FUNCTION);
            break;
          }
        }
      }
      else {

        /*--- Build the grid data structures using the ParMETIS coloring. ---*/
        
        geometry_container[iZone][iInst][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);

      }

      /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

      delete geometry_aux;

      /*--- Add the Send/Receive boundaries ---*/
      geometry_container[iZone][iInst][MESH_0]->SetSendReceive(config_container[iZone]);

      /*--- Add the Send/Receive boundaries ---*/
      geometry_container[iZone][iInst][MESH_0]->SetBoundaries(config_container[iZone]);

    }

  }

}

void CErrorEstimationDriver::Geometrical_Preprocessing(CConfig **config_container, CGeometry ****geometry_container) {

  bool fea = false;

  for (iZone = 0; iZone < nZone; iZone++) {

    fea = ((config_container[iZone]->GetKind_Solver() == FEM_ELASTICITY) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM));

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Compute elements surrounding points, points surrounding points ---*/

      if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetPoint_Connectivity();

      /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/

      // if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
      // geometry_container[iZone][iInst][MESH_0]->SetRCM_Ordering(config_container[iZone]);

      /*--- recompute elements surrounding points, points surrounding points ---*/

      if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetPoint_Connectivity();

      /*--- Compute elements surrounding elements ---*/

      if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetElement_Connectivity();

      /*--- Check the orientation before computing geometrical quantities ---*/

      geometry_container[iZone][iInst][MESH_0]->SetBoundVolume();
      if (config_container[iZone]->GetReorientElements()) {
        if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
        geometry_container[iZone][iInst][MESH_0]->Check_IntElem_Orientation(config_container[iZone]);
        geometry_container[iZone][iInst][MESH_0]->Check_BoundElem_Orientation(config_container[iZone]);
      }

      /*--- Create the edge structure ---*/

      if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetEdges();
      geometry_container[iZone][iInst][MESH_0]->SetVertex(config_container[iZone]);

      /*--- Compute cell center of gravity ---*/

      if ((rank == MASTER_NODE) && (!fea)) cout << "Computing centers of gravity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetCoord_CG();

      /*--- Create the control volume structures ---*/

      if ((rank == MASTER_NODE) && (!fea)) cout << "Setting the control volume structure." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetControlVolume(config_container[iZone], ALLOCATE);
      geometry_container[iZone][iInst][MESH_0]->SetBoundControlVolume(config_container[iZone], ALLOCATE);

      /*--- Visualize a dual control volume if requested ---*/

      if ((config_container[iZone]->GetVisualize_CV() >= 0) &&
          (config_container[iZone]->GetVisualize_CV() < (long)geometry_container[iZone][iInst][MESH_0]->GetnPointDomain()))
        geometry_container[iZone][iInst][MESH_0]->VisualizeControlVolume(config_container[iZone], UPDATE);

      /*--- Identify closest normal neighbor ---*/

      if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
      geometry_container[iZone][iInst][MESH_0]->FindNormal_Neighbor(config_container[iZone]);

      /*--- Store the global to local mapping. ---*/

      if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetGlobal_to_Local_Point();

      /*--- Compute the surface curvature ---*/

      if ((rank == MASTER_NODE) && (!fea)) cout << "Compute the surface curvature." << endl;
      geometry_container[iZone][iInst][MESH_0]->ComputeSurf_Curvature(config_container[iZone]);

      /*--- Check for periodicity and disable MG if necessary. ---*/

      if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
      geometry_container[iZone][iInst][MESH_0]->Check_Periodicity(config_container[iZone]);

    }

  }

  /*--- Create the data structure for MPI point-to-point communications. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        geometry_container[iZone][iInst][MESH_0]->PreprocessP2PComms(geometry_container[iZone][iInst][MESH_0], config_container[iZone]);
    }
  }
  
  /*--- Perform a few preprocessing routines and communications. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        
      /*--- Compute the max length. ---*/
        
      if ((rank == MASTER_NODE) && (!fea)) cout << "Finding max control volume width." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetMaxLength(config_container[iZone]);
        
      /*--- Communicate the number of neighbors. This is needed for
       some centered schemes and for multigrid in parallel. ---*/
        
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (!fea)) cout << "Communicating number of neighbors." << endl;
      geometry_container[iZone][iInst][MESH_0]->InitiateComms(geometry_container[iZone][iInst][MESH_0], config_container[iZone], NEIGHBORS);
      geometry_container[iZone][iInst][MESH_0]->CompleteComms(geometry_container[iZone][iInst][MESH_0], config_container[iZone], NEIGHBORS);
    }
  }

}

void CErrorEstimationDriver::Geometrical_Preprocessing_DGFEM(CConfig **config_container, CGeometry ****geometry_container) {

  //*--- Loop over the number of zones of the fine grid. ---*/

  for(unsigned short iZone = 0; iZone < nZone; iZone++) {

    /*--- Loop over the time instances of this zone. ---*/
    for(unsigned short iInst = 0; iInst < nInst[iZone]; iInst++) {

      /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
       define all virtual functions in the base class CGeometry. ---*/
      CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry_container[iZone][iInst][MESH_0]);

      /*--- Determine the standard elements for the volume elements. ---*/
      if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
      DGMesh->CreateStandardVolumeElements(config_container[iZone]);

      /*--- Create the face information needed to compute the contour integral
       for the elements in the Discontinuous Galerkin formulation. ---*/
      if (rank == MASTER_NODE) cout << "Creating face information." << endl;
      DGMesh->CreateFaces(config_container[iZone]);

      /*--- Compute the metric terms of the volume elements. ---*/
      if (rank == MASTER_NODE) cout << "Computing metric terms volume elements." << endl;
      DGMesh->MetricTermsVolumeElements(config_container[iZone]);

      /*--- Compute the metric terms of the surface elements. ---*/
      if (rank == MASTER_NODE) cout << "Computing metric terms surface elements." << endl;
      DGMesh->MetricTermsSurfaceElements(config_container[iZone]);

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

      /*--- Initialize the static mesh movement, if necessary. ---*/
      const unsigned short Kind_Grid_Movement = config_container[iZone]->GetKind_GridMovement(iZone);
      const bool initStaticMovement = (config_container[iZone]->GetGrid_Movement() &&
                                      (Kind_Grid_Movement == MOVING_WALL    ||
                                       Kind_Grid_Movement == ROTATING_FRAME ||
                                       Kind_Grid_Movement == STEADY_TRANSLATION));

      if(initStaticMovement){
        if (rank == MASTER_NODE) cout << "Initialize Static Mesh Movement" << endl;
        DGMesh->InitStaticMeshMovement(config_container[iZone], Kind_Grid_Movement, iZone);
      }

      /*--- Perform the preprocessing tasks when wall functions are used. ---*/
      if (rank == MASTER_NODE) cout << "Preprocessing for the wall functions. " << endl;
      DGMesh->WallFunctionPreprocessing(config_container[iZone]);

      /*--- Store the global to local mapping. ---*/
      if (rank == MASTER_NODE) cout << "Storing a mapping from global to local DOF index." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetGlobal_to_Local_Point();
    }

    /*--- Loop to create the coarser grid levels. ---*/

    for(unsigned short iMGlevel=1; iMGlevel<=config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

      SU2_MPI::Error("Geometrical_Preprocessing_DGFEM: Coarse grid levels not implemented yet.",
                     CURRENT_FUNCTION);
    }
  }
}

void CErrorEstimationDriver::Solver_Preprocessing(CSolver ****solver_container, CGeometry ***geometry,
                                                   CConfig *config, unsigned short val_iInst) {
  
  bool euler, ns, turbulent,
  fem_euler, fem_ns, fem_turbulent, fem_transition,
  adj_euler, adj_ns, adj_turb,
  heat_fvm,
  fem, disc_adj_fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_heat,
  fem_dg_flow, fem_dg_shock_persson,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;  

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent     = false;
  fem_euler        = false;  fem_ns          = false;  fem_turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb      = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat_fvm         = false;  disc_adj_heat    = false;
  transition       = false;  fem_transition   = false;
  template_solver  = false;
  fem_dg_flow      = false;  fem_dg_shock_persson = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;
  
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case FEM_LES : fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Determine the kind of FEM solver used for the flow. ---*/

  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }

  /*--- Determine the kind of shock capturing method for FEM DG solver. ---*/

  switch( config->GetKind_FEM_DG_Shock() ) {
    case PERSSON: fem_dg_shock_persson = true; break;
  }

  /*--- Assign turbulence model booleans ---*/

  if (turbulent || fem_turbulent){
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;           break;
      case SA_E:   e_spalart_allmaras = true;   break;
      case SA_COMP: comp_spalart_allmaras = true; break;
      case SA_E_COMP: e_comp_spalart_allmaras = true; break;
      default: SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION); break;
    }
  }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][INSTANCE][MESH_0][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/


    /*--- Allocate solution for a template problem ---*/

  if (template_solver) {
    solver_container[val_iInst][MESH_0][TEMPLATE_SOL] = new CTemplateSolver(geometry[val_iInst][MESH_0], config);
    DOFsPerPoint += solver_container[val_iInst][MESH_0][TEMPLATE_SOL]->GetnVar();
  }

    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
    /*--- Note that we need both direct and adjoint for the error estimate ---*/

  if (euler || disc_adj) {
    if (compressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CEulerSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    if (incompressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CIncEulerSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    DOFsPerPoint += solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  }
  if (ns) {
    if (compressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CNSSolver(geometry[val_iInst][MESH_0], config, MESH_0);
    }
    if (incompressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CIncNSSolver(geometry[val_iInst][MESH_0], config, MESH_0);
    }
    DOFsPerPoint += solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  }
  if (turbulent) {
    if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras || neg_spalart_allmaras) {
      solver_container[val_iInst][MESH_0][TURB_SOL] = new CTurbSASolver(geometry[val_iInst][MESH_0], config, MESH_0, solver_container[val_iInst][MESH_0][FLOW_SOL]->GetFluidModel() );
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      solver_container[val_iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0);
    }
    else if (menter_sst) {
      solver_container[val_iInst][MESH_0][TURB_SOL] = new CTurbSSTSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      solver_container[val_iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    DOFsPerPoint += solver_container[val_iInst][MESH_0][TURB_SOL]->GetnVar();
    if (transition) {
      solver_container[val_iInst][MESH_0][TRANS_SOL] = new CTransLMSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      DOFsPerPoint += solver_container[val_iInst][MESH_0][TRANS_SOL]->GetnVar();
    }
  }

    /*--- Allocate solution for adjoint problem ---*/

  if (disc_adj || euler || ns || turbulent) {
    solver_container[val_iInst][MESH_0][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[val_iInst][MESH_0], config, solver_container[val_iInst][MESH_0][FLOW_SOL], RUNTIME_FLOW_SYS, MESH_0);
    DOFsPerPoint += solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->GetnVar();
    if (disc_adj_turb || turbulent) {
      solver_container[val_iInst][MESH_0][ADJTURB_SOL] = new CDiscAdjSolver(geometry[val_iInst][MESH_0], config, solver_container[val_iInst][MESH_0][TURB_SOL], RUNTIME_TURB_SYS, MESH_0);
      DOFsPerPoint += solver_container[val_iInst][MESH_0][ADJTURB_SOL]->GetnVar();
    }
  }  

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool update_geo = true;
  if (config->GetFSI_Simulation()) update_geo = false;

  Solver_Restart(solver_container, geometry, config, update_geo, val_iInst);

}

void CErrorEstimationDriver::Solver_Restart(CSolver ****solver_container, CGeometry ***geometry,
                             CConfig *config, bool update_geo, unsigned short val_iInst) {

  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat_fvm, fem, fem_euler, fem_ns, fem_dg_flow,
  template_solver, disc_adj, disc_adj_fem, disc_adj_turb, disc_adj_heat;
  int val_iter = 0;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns           = false;  turbulent   = false;
  adj_euler        = false;  adj_ns       = false;  adj_turb    = false;
  fem_euler        = false;  fem_ns       = false;  fem_dg_flow = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  disc_adj_turb    = false;
  heat_fvm         = false;  disc_adj_heat    = false;
  template_solver  = false;

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool restart      = config->GetRestart();
  bool restart_flow = config->GetRestart_Flow();
  bool no_restart   = false;

  /*--- Adjust iteration number for unsteady restarts. ---*/

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC); // Dynamic simulation (FSI).

  if (dual_time) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
  }

  if (time_stepping) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
  }

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; break;
    case FEM_LES : fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;

  }

  /*--- Determine the kind of FEM solver used for the flow. ---*/

  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }

  /*--- Load restarts for any of the active solver containers. Note that
   these restart routines fill the fine grid and interpolate to all MG levels. ---*/

  if (restart || restart_flow) {
    if (euler || ns) {
      solver_container[val_iInst][MESH_0][FLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (turbulent) {
      solver_container[val_iInst][MESH_0][TURB_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (fem) {
      if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
      solver_container[val_iInst][MESH_0][FEA_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (fem_euler || fem_ns) {
      if (fem_dg_flow)
        solver_container[val_iInst][MESH_0][FLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (heat_fvm) {
      solver_container[val_iInst][MESH_0][HEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
  }

  if (restart) {
    if (template_solver) {
      no_restart = true;
    }
    if (heat_fvm) {
      solver_container[val_iInst][MESH_0][HEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (adj_euler || adj_ns) {
      solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (adj_turb) {
      no_restart = true;
    }
    if (disc_adj) {
      solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
      if (disc_adj_turb)
        solver_container[val_iInst][MESH_0][ADJTURB_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
      if (disc_adj_heat)
        solver_container[val_iInst][MESH_0][ADJHEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (disc_adj_fem) {
        if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
        solver_container[val_iInst][MESH_0][ADJFEA_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (disc_adj_heat) {
      solver_container[val_iInst][MESH_0][ADJHEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
  }

  /*--- Exit if a restart was requested for a solver that is not available. ---*/

  if (no_restart) {
    SU2_MPI::Error(string("A restart capability has not been implemented yet for this solver.\n") +
                   string("Please set RESTART_SOL= NO and try again."), CURRENT_FUNCTION);
  }

  /*--- Think about calls to pre / post-processing here, plus realizability checks. ---*/
  
}

void CErrorEstimationDriver::ComputeMetric() {

  CSolver *solver_flow = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL],
          *solver_adj  = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];

  if (rank == MASTER_NODE)
    cout << endl <<"----------------------------- Compute Metric ----------------------------" << endl;

  //--- 2D
  if(nDim == 2){
    //--- Volume flow grad
    if(rank == MASTER_NODE) cout << "Computing flow volume gradient via L2 Projection." << endl;
    solver_flow->SetGradient_L2Proj2(geometry_container[ZONE_0][INST_0][MESH_0], 
                                     config_container[ZONE_0]);

    //--- Volume flow Hess
    if(rank == MASTER_NODE) cout << "Computing flow volume Hessian via L2 Projection." << endl;
    solver_flow->SetHessian_L2Proj2(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);

    //--- Volume adj grad
    if(rank == MASTER_NODE) cout << "Computing adjoint volume gradient via L2 Projection." << endl;
    solver_adj->SetGradient_L2Proj2(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);

    //--- Metric
    if(rank == MASTER_NODE) cout << "Computing goal-oriented metric tensor." << endl;
    SumWeightedHessian2(solver_flow, solver_adj, geometry_container[ZONE_0][INST_0][MESH_0]);
  }

  //--- 3D
  else{
    //--- Volume flow grad
    if(rank == MASTER_NODE) cout << "Computing flow volume gradient via L2 Projection." << endl;
    solver_flow->SetGradient_L2Proj3(geometry_container[ZONE_0][INST_0][MESH_0], 
                                     config_container[ZONE_0]);

    //--- Volume flow Hess
    if(rank == MASTER_NODE) cout << "Computing flow volume Hessian via L2 Projection." << endl;
    solver_flow->SetHessian_L2Proj3(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);

    //--- Volume adj grad
    if(rank == MASTER_NODE) cout << "Computing adjoint volume gradient via L2 Projection." << endl;
    solver_adj->SetGradient_L2Proj3(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);

    //--- Metric
    if(rank == MASTER_NODE) cout << "Computing goal-oriented metric tensor." << endl;
    SumWeightedHessian3(solver_flow, solver_adj, geometry_container[ZONE_0][INST_0][MESH_0]);
  }
}

void CErrorEstimationDriver::SumWeightedHessian2(CSolver   *solver_flow,
                                                 CSolver   *solver_adj,
                                                 CGeometry *geometry) {

  unsigned long iPoint, nPointDomain = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain();
  unsigned short nVarMetr = solver_flow->GetnVar(), 
                 nFluxMetr = 2;  //--- TODO: adjust size for goal (for different solvers, currently Euler) vs. feature
  unsigned short nMetr = 3*(nDim-1);

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    //--- initialize metric tensor
    for(unsigned short im = 0; im < nMetr; ++im)
      solver_flow->node[iPoint]->SetAnisoMetr(im, 0.0);

    //--- perform summation of weighted Hessians
    for (unsigned short iVar = 0; iVar < nVarMetr; ++iVar) {

      for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
        const unsigned short ig = iVar*nDim + iFlux;
        const su2double grad = solver_adj->node[iPoint]->GetAnisoGrad(ig);

        for (unsigned short im = 0; im < nMetr; ++im) {
          const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
          const su2double hess = solver_flow->node[iPoint]->GetAnisoHess(ih);
          const su2double part = abs(grad)*hess;
          solver_flow->node[iPoint]->AddAnisoMetr(im,part);
        }
      }
    }
  }

  //--- avoid null metrics and obtain global complexity
  su2double localComplexity  = 0.0,
            globalComplexity = 0.0,
            p = 1.0,                                                                    // For now, hardcode L1 metric
            outComplexity = su2double(config_container[ZONE_0]->GetMesh_Complexity());  // Constraint mesh complexity
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    CVariable *var = solver_flow->node[iPoint];

    const su2double a = var->GetAnisoMetr(0);
    const su2double b = var->GetAnisoMetr(1);
    const su2double c = var->GetAnisoMetr(2);
    const su2double d = pow(a-c,2.) + 4.*b*b;

    const su2double Lam1 = (a+c+sqrt(d))/2.;
    const su2double Lam2 = (a+c-sqrt(d))/2.;

    su2double RuH[2][2] = {{b, b},
                           {Lam1-a, Lam2-a}};

    if(abs(b) < 1.0e-16){
      RuH[0][0] = 1.0; RuH[0][1] = 0.0;
      RuH[1][0] = 0.0; RuH[1][1] = 1.0;
    }

    const su2double RuU[2][2]    = {{RuH[0][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]), RuH[0][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1])},
                                    {RuH[1][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]), RuH[1][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1])}};

    const su2double Lam1new = max(abs(Lam1), 1.0E-16);
    const su2double Lam2new = max(abs(Lam2), 1.0E-16);

    const su2double LamRuU[2][2] = {{abs(Lam1new)*RuU[0][0],abs(Lam1new)*RuU[1][0]},
                                    {abs(Lam2new)*RuU[0][1],abs(Lam2new)*RuU[1][1]}};

    const su2double MetrNew[3]   = {RuU[0][0]*LamRuU[0][0]+RuU[0][1]*LamRuU[1][0], 
                                    RuU[0][0]*LamRuU[0][1]+RuU[0][1]*LamRuU[1][1], 
                                    RuU[1][0]*LamRuU[0][1]+RuU[1][1]*LamRuU[1][1]};

    var->SetAnisoMetr(0, MetrNew[0]);
    var->SetAnisoMetr(1, MetrNew[1]);
    var->SetAnisoMetr(2, MetrNew[2]);

    const su2double Vol = geometry->node[iPoint]->GetVolume();

    localComplexity += pow(abs(Lam1*Lam2),p/(2.*p+3.))*Vol;
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&localComplexity, &globalComplexity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  globalComplexity = localComplexity;
#endif

  if(rank == MASTER_NODE)
    cout << "Current mesh complexity: " << globalComplexity << "." << endl;

  //--- normalize to obtain Lp metric then constrain size
  su2double hmax = config_container[ZONE_0]->GetMesh_Hmax(),
            hmin = config_container[ZONE_0]->GetMesh_Hmin();
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    CVariable *var = solver_flow->node[iPoint];

    const su2double a = var->GetAnisoMetr(0);
    const su2double b = var->GetAnisoMetr(1);
    const su2double c = var->GetAnisoMetr(2);
    const su2double d = pow(a-c,2.) + 4.*b*b;

    const su2double Lam1 = (a+c+sqrt(d))/2.;
    const su2double Lam2 = (a+c-sqrt(d))/2.;

    const su2double factor = pow(outComplexity/globalComplexity, 2./3.)
                           * pow(abs(Lam1*Lam2), -1./(2.*p+3.));

    su2double RuH[2][2] = {{b, b},
                           {Lam1-a, Lam2-a}};

    if(abs(b) < 1.0e-16){
      RuH[0][0] = 1.0; RuH[0][1] = 0.0;
      RuH[1][0] = 0.0; RuH[1][1] = 1.0;
    }

    const su2double RuU[2][2]    = {{RuH[0][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]), RuH[0][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1])},
                                    {RuH[1][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]), RuH[1][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1])}};

    const su2double Lam1new = min(max(abs(factor*Lam1), 1./(hmax*hmax)), 1./(hmin*hmin));
    const su2double Lam2new = min(max(abs(factor*Lam2), 1./(hmax*hmax)), 1./(hmin*hmin));

    const su2double LamRuU[2][2] = {{abs(Lam1new)*RuU[0][0],abs(Lam1new)*RuU[1][0]},
                                    {abs(Lam2new)*RuU[0][1],abs(Lam2new)*RuU[1][1]}};

    const su2double MetrNew[3]   = {RuU[0][0]*LamRuU[0][0]+RuU[0][1]*LamRuU[1][0], 
                                    RuU[0][0]*LamRuU[0][1]+RuU[0][1]*LamRuU[1][1], 
                                    RuU[1][0]*LamRuU[0][1]+RuU[1][1]*LamRuU[1][1]};

    var->SetAnisoMetr(0, MetrNew[0]);
    var->SetAnisoMetr(1, MetrNew[1]);
    var->SetAnisoMetr(2, MetrNew[2]);
  }
}

void CErrorEstimationDriver::SumWeightedHessian3(CSolver   *solver_flow,
                                                 CSolver   *solver_adj,
                                                 CGeometry *geometry) {

  unsigned long iPoint, nPointDomain = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain();
  unsigned short nVarMetr = solver_flow->GetnVar(), 
                 nFluxMetr = 3;  //--- TODO: adjust size for goal (for different solvers, currently Euler) vs. feature
  unsigned short nMetr = 3*(nDim-1);

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    //--- initialize metric tensor
    for(unsigned short im = 0; im < nMetr; ++im)
      solver_flow->node[iPoint]->SetAnisoMetr(im, 0.0);

    //--- perform summation of weighted Hessians
    for (unsigned short iVar = 0; iVar < nVarMetr; ++iVar) {

      for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
        const unsigned short ig = iVar*nDim + iFlux;
        const su2double grad = solver_adj->node[iPoint]->GetAnisoGrad(ig);

        for (unsigned short im = 0; im < nMetr; ++im) {
          const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
          const su2double hess = solver_flow->node[iPoint]->GetAnisoHess(ih);
          const su2double part = abs(grad)*hess;
          solver_flow->node[iPoint]->AddAnisoMetr(im,part);
        }
      }
    }
  }

  //--- avoid null metrics and obtain global complexity
  su2double localComplexity  = 0.0,
            globalComplexity = 0.0,
            p = 1.0,                                                                    // For now, hardcode L1 metric
            outComplexity = su2double(config_container[ZONE_0]->GetMesh_Complexity());  // Constraint mesh complexity
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    CVariable *var = solver_flow->node[iPoint];

    su2double Lam[3], RuH[3][3];

    const su2double a = var->GetAnisoMetr(0);
    const su2double b = var->GetAnisoMetr(1);
    const su2double c = var->GetAnisoMetr(2);
    const su2double d = var->GetAnisoMetr(3);
    const su2double e = var->GetAnisoMetr(4);
    const su2double f = var->GetAnisoMetr(5);

    const su2double p1 = b*b + c*c + e*e;
    if(p1 < 1.0e-16){
      //--- eigenvalues
      Lam[0] = abs(a);
      Lam[1] = abs(d);
      Lam[2] = abs(f);

      //--- eigenvectors
      RuH[0][0] = RuH[1][1] = RuH[2][2] = 1.;
      RuH[0][1] = RuH[1][0] = RuH[0][2] = RuH[2][0] = RuH[1][2] = RuH[2][1] = 0.;
    }
    else{
      //--- eigenvalues
      const su2double q  = (a + d + f)/3.;
      const su2double p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2.*p1;
      const su2double p3 = sqrt(p2/6.);

      const su2double aa = (a-q)/p3, bb = b/p3, cc = c/p3,
                      dd = (d-q)/p3, ee = e/p3, ff = (f-q)/p3;

      const su2double det = 0.5*(aa*(dd*ff - ee*ee) + bb*(cc*ee - bb*ff) + cc*(bb*ee - dd*cc));

      const su2double pi = 3.141592653589793238;
      su2double phi;
      if(det <= -1.)     phi = pi/3.;
      else if(det >= 1.) phi = 0.;
      else               phi = acos(det)/3.;

      Lam[0] = q+2.*p*cos(phi);
      Lam[1] = q+2.*p*cos(phi+2.*p3/3.);
      Lam[2] = 3.*q-Lam[0]-Lam[1];

      //--- eigenvectors
      for(unsigned short i = 0; i < 3; ++i) {
        su2double A[3][4] = {{a-Lam[i], b,        c,        0.},
                             {b,        d-Lam[i], e,        0.},
                             {c,        e,        f-Lam[i], 0.}};

        int nr = 3, nc = 4, j = 0;
        while(j < nr){
          for(int r = 0; r < nr; ++r) {
            const su2double d = A[j][j];
            const su2double m = A[r][j]/d;
            for(int c = 0; c < nc; ++c) {
              if(r == j) A[r][c] /= d;
              else       A[r][c] -= A[j][c]*m;
            }
          }
          j++;
        }
        RuH[0][i] = A[0][4];
        RuH[1][i] = A[1][4];
        RuH[2][i] = A[2][4];
      }

      Lam[0] = abs(Lam[0]);
      Lam[1] = abs(Lam[1]);
      Lam[2] = abs(Lam[2]);
    }

    const su2double RuU[3][3]    = {{RuH[0][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]+RuH[2][0]*RuH[2][0]),
                                     RuH[0][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1]+RuH[2][1]*RuH[2][1]),
                                     RuH[0][2]/sqrt(RuH[0][2]*RuH[0][2]+RuH[1][2]*RuH[1][2]+RuH[2][2]*RuH[2][2])},
                                    {RuH[1][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]+RuH[2][0]*RuH[2][0]),
                                     RuH[1][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1]+RuH[2][1]*RuH[2][1]),
                                     RuH[1][2]/sqrt(RuH[0][2]*RuH[0][2]+RuH[1][2]*RuH[1][2]+RuH[2][2]*RuH[2][2])},
                                    {RuH[2][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]+RuH[2][0]*RuH[2][0]),
                                     RuH[2][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1]+RuH[2][1]*RuH[2][1]),
                                     RuH[2][2]/sqrt(RuH[0][2]*RuH[0][2]+RuH[1][2]*RuH[1][2]+RuH[2][2]*RuH[2][2])}};

    const su2double LamRuU[3][3] = {{Lam[0]*RuU[0][0],Lam[0]*RuU[1][0],Lam[0]*RuU[2][0]},
                                    {Lam[1]*RuU[0][1],Lam[1]*RuU[1][1],Lam[1]*RuU[2][1]},
                                    {Lam[2]*RuU[0][2],Lam[2]*RuU[1][2],Lam[2]*RuU[2][2]}};

    const su2double Metr[6]      = {RuU[0][0]*LamRuU[0][0]+RuU[0][1]*LamRuU[1][0]+RuU[0][2]*LamRuU[2][0], 
                                    RuU[0][0]*LamRuU[0][1]+RuU[0][1]*LamRuU[1][1]+RuU[0][2]*LamRuU[2][1],
                                    RuU[0][0]*LamRuU[0][2]+RuU[0][1]*LamRuU[1][2]+RuU[0][2]*LamRuU[2][2],
                                    RuU[1][0]*LamRuU[0][1]+RuU[1][1]*LamRuU[1][1]+RuU[1][2]*LamRuU[2][1],
                                    RuU[1][0]*LamRuU[0][2]+RuU[1][1]*LamRuU[1][2]+RuU[1][2]*LamRuU[2][2],
                                    RuU[2][0]*LamRuU[0][2]+RuU[2][1]*LamRuU[1][2]+RuU[2][2]*LamRuU[2][2]};

    var->SetAnisoMetr(0, Metr[0]);
    var->SetAnisoMetr(1, Metr[1]);
    var->SetAnisoMetr(2, Metr[2]);
    var->SetAnisoMetr(3, Metr[3]);
    var->SetAnisoMetr(4, Metr[4]);
    var->SetAnisoMetr(5, Metr[5]);

    const su2double Vol = geometry->node[iPoint]->GetVolume();

    localComplexity += pow(abs(Lam[0]*Lam[1]*Lam[2]),p/(2.*p+3.))*Vol;
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&localComplexity, &globalComplexity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  globalComplexity = localComplexity;
#endif

  if(rank == MASTER_NODE)
    cout << "Current mesh complexity: " << globalComplexity << "." << endl;

  //--- normalize to obtain Lp metric then constrain size
  su2double hmax = config_container[ZONE_0]->GetMesh_Hmax(),
            hmin = config_container[ZONE_0]->GetMesh_Hmin();
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    CVariable *var = solver_flow->node[iPoint];

    su2double Lam[3], RuH[3][3];

    const su2double a = var->GetAnisoMetr(0);
    const su2double b = var->GetAnisoMetr(1);
    const su2double c = var->GetAnisoMetr(2);
    const su2double d = var->GetAnisoMetr(3);
    const su2double e = var->GetAnisoMetr(4);
    const su2double f = var->GetAnisoMetr(5);

    const su2double p1 = b*b + c*c + e*e;
    if(p1 < 1.0e-16){
      //--- eigenvalues
      Lam[0] = abs(a);
      Lam[1] = abs(d);
      Lam[2] = abs(f);

      //--- eigenvectors
      RuH[0][0] = RuH[1][1] = RuH[2][2] = 1.;
      RuH[0][1] = RuH[1][0] = RuH[0][2] = RuH[2][0] = RuH[1][2] = RuH[2][1] = 0.;
    }
    else{
      //--- eigenvalues
      const su2double q  = (a + d + f)/3.;
      const su2double p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2.*p1;
      const su2double p3 = sqrt(p2/6.);

      const su2double aa = (a-q)/p3, bb = b/p3, cc = c/p3,
                      dd = (d-q)/p3, ee = e/p3, ff = (f-q)/p3;

      const su2double det = 0.5*(aa*(dd*ff - ee*ee) + bb*(cc*ee - bb*ff) + cc*(bb*ee - dd*cc));

      const su2double pi = 3.141592653589793238;
      su2double phi;
      if(det <= -1.)     phi = pi/3.;
      else if(det >= 1.) phi = 0.;
      else               phi = acos(det)/3.;

      Lam[0] = q+2.*p*cos(phi);
      Lam[1] = q+2.*p*cos(phi+2.*p3/3.);
      Lam[2] = 3.*q-Lam[0]-Lam[1];

      //--- eigenvectors
      for(unsigned short i = 0; i < 3; ++i) {
        su2double A[3][4] = {{a-Lam[i], b,        c,        0.},
                             {b,        d-Lam[i], e,        0.},
                             {c,        e,        f-Lam[i], 0.}};

        int nr = 3, nc = 4, j = 0;
        while(j < nr){
          for(int r = 0; r < nr; ++r) {
            const su2double d = A[j][j];
            const su2double m = A[r][j]/d;
            for(int c = 0; c < nc; ++c) {
              if(r == j) A[r][c] /= d;
              else       A[r][c] -= A[j][c]*m;
            }
          }
          j++;
        }
        RuH[0][i] = A[0][4];
        RuH[1][i] = A[1][4];
        RuH[2][i] = A[2][4];
      }

      Lam[0] = abs(Lam[0]);
      Lam[1] = abs(Lam[1]);
      Lam[2] = abs(Lam[2]);
    }

    const su2double factor = pow(outComplexity/globalComplexity, 2./3.)
                           * pow(abs(Lam[0]*Lam[1]*Lam[2]), -1./(2.*p+3.));

    const su2double Lam1new = min(max(abs(factor*Lam[0]), 1./(hmax*hmax)), 1./(hmin*hmin));
    const su2double Lam2new = min(max(abs(factor*Lam[1]), 1./(hmax*hmax)), 1./(hmin*hmin));
    const su2double Lam3new = min(max(abs(factor*Lam[2]), 1./(hmax*hmax)), 1./(hmin*hmin));

    const su2double RuU[3][3]    = {{RuH[0][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]+RuH[2][0]*RuH[2][0]),
                                     RuH[0][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1]+RuH[2][1]*RuH[2][1]),
                                     RuH[0][2]/sqrt(RuH[0][2]*RuH[0][2]+RuH[1][2]*RuH[1][2]+RuH[2][2]*RuH[2][2])},
                                    {RuH[1][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]+RuH[2][0]*RuH[2][0]),
                                     RuH[1][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1]+RuH[2][1]*RuH[2][1]),
                                     RuH[1][2]/sqrt(RuH[0][2]*RuH[0][2]+RuH[1][2]*RuH[1][2]+RuH[2][2]*RuH[2][2])},
                                    {RuH[2][0]/sqrt(RuH[0][0]*RuH[0][0]+RuH[1][0]*RuH[1][0]+RuH[2][0]*RuH[2][0]),
                                     RuH[2][1]/sqrt(RuH[0][1]*RuH[0][1]+RuH[1][1]*RuH[1][1]+RuH[2][1]*RuH[2][1]),
                                     RuH[2][2]/sqrt(RuH[0][2]*RuH[0][2]+RuH[1][2]*RuH[1][2]+RuH[2][2]*RuH[2][2])}};

    const su2double LamRuU[3][3] = {{Lam1new*RuU[0][0],Lam1new*RuU[1][0],Lam1new*RuU[2][0]},
                                    {Lam2new*RuU[0][1],Lam2new*RuU[1][1],Lam2new*RuU[2][1]},
                                    {Lam3new*RuU[0][2],Lam3new*RuU[1][2],Lam3new*RuU[2][2]}};

    const su2double MetrNew[6]   = {RuU[0][0]*LamRuU[0][0]+RuU[0][1]*LamRuU[1][0]+RuU[0][2]*LamRuU[2][0], 
                                    RuU[0][0]*LamRuU[0][1]+RuU[0][1]*LamRuU[1][1]+RuU[0][2]*LamRuU[2][1],
                                    RuU[0][0]*LamRuU[0][2]+RuU[0][1]*LamRuU[1][2]+RuU[0][2]*LamRuU[2][2],
                                    RuU[1][0]*LamRuU[0][1]+RuU[1][1]*LamRuU[1][1]+RuU[1][2]*LamRuU[2][1],
                                    RuU[1][0]*LamRuU[0][2]+RuU[1][1]*LamRuU[1][2]+RuU[1][2]*LamRuU[2][2],
                                    RuU[2][0]*LamRuU[0][2]+RuU[2][1]*LamRuU[1][2]+RuU[2][2]*LamRuU[2][2]};

    //--- swap indices 2 and 3 to be consistent with AMG, which stores lower triangle in 3D
    var->SetAnisoMetr(0, MetrNew[0]);
    var->SetAnisoMetr(1, MetrNew[1]);
    var->SetAnisoMetr(2, MetrNew[3]);
    var->SetAnisoMetr(3, MetrNew[2]);
    var->SetAnisoMetr(4, MetrNew[4]);
    var->SetAnisoMetr(5, MetrNew[5]);
  }
}

void CErrorEstimationDriver::Output() {

  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter_OffSet();

  /*--- Set Kind_Solver to primal since adaptation requires Mach (or pressure) ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
    case DISC_ADJ_EULER: config_container[ZONE_0]->SetKind_Solver(EULER); break;
    case DISC_ADJ_NAVIER_STOKES: config_container[ZONE_0]->SetKind_Solver(NAVIER_STOKES); break;
    case DISC_ADJ_RANS: config_container[ZONE_0]->SetKind_Solver(RANS); break;
  }

  /*--- Set DiscreteAdjoint flag to false so we write to correct files ---*/

  config_container[ZONE_0]->SetDiscrete_Adjoint(false);

  if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------";

    /*--- Execute the routine for writing restart, volume solution,
     surface solution, and surface comma-separated value files. ---*/

  output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, ExtIter, nZone);

  if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;

  
}

void CErrorEstimationDriver::Postprocessing() {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Error Estimation Postprocessing ---------------------" << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      Solver_Postprocessing(solver_container[iZone], geometry_container[iZone][iInst],
                            config_container[iZone], iInst);
    } // End of loop over iInst
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver containers." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone] != NULL) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        for (unsigned short iMGlevel = 0; iMGlevel < config_container[iZone]->GetnMGLevels()+1; iMGlevel++) {
          if (geometry_container[iZone][iInst][iMGlevel] != NULL) delete geometry_container[iZone][iInst][iMGlevel];
        }
        if (geometry_container[iZone][iInst] != NULL) delete [] geometry_container[iZone][iInst];
      } // End of loop over iInst
      delete [] geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  if (rank == MASTER_NODE) cout << "Deleted CGeometry containers." << endl;

  /*--- Deallocate config container ---*/
  if (config_container!= NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig containers." << endl;

}

void CErrorEstimationDriver::Solver_Postprocessing(CSolver ****solver_container, CGeometry **geometry,
                                                   CConfig *config, unsigned short val_iInst) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat_fvm, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_fem, disc_adj_heat,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj        = false;
  fem              = false;  disc_adj_fem    = false;
  heat_fvm        = false;   disc_adj_heat   = false;
  transition       = false;
  template_solver  = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : euler = true; break;
    case FEM_NAVIER_STOKES:
    case FEM_LES: ns = true; break;
    case FEM_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
    case SA:     spalart_allmaras = true;     break;
    case SA_NEG: neg_spalart_allmaras = true; break;
    case SST:    menter_sst = true;           break;
    case SA_E: e_spalart_allmaras = true; break;
    case SA_COMP: comp_spalart_allmaras = true; break;
    case SA_E_COMP: e_comp_spalart_allmaras = true; break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    
    /*--- DeAllocate solution for a template problem ---*/
    
    if (template_solver) {
      delete solver_container[val_iInst][iMGlevel][TEMPLATE_SOL];
    }

    /*--- DeAllocate solution for adjoint problem ---*/
    
    if (adj_euler || adj_ns || disc_adj) {
      delete solver_container[val_iInst][iMGlevel][ADJFLOW_SOL];
      if (disc_adj_turb || adj_turb) {
        delete solver_container[val_iInst][iMGlevel][ADJTURB_SOL];
      }
      if (heat_fvm) {
        delete solver_container[val_iInst][iMGlevel][ADJHEAT_SOL];
      }
    }

    if (disc_adj_heat) {
      delete solver_container[val_iInst][iMGlevel][ADJHEAT_SOL];
    }

    /*--- DeAllocate solution for direct problem ---*/
    
    if (euler || ns) {
      delete solver_container[val_iInst][iMGlevel][FLOW_SOL];
    }

    if (turbulent) {
      if (spalart_allmaras || neg_spalart_allmaras || menter_sst || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
        delete solver_container[val_iInst][iMGlevel][TURB_SOL];
      }
      if (transition) {
        delete solver_container[val_iInst][iMGlevel][TRANS_SOL];
      }
    }
    if (heat_fvm) {
      delete solver_container[val_iInst][iMGlevel][HEAT_SOL];
    }
    if (fem) {
      delete solver_container[val_iInst][iMGlevel][FEA_SOL];
    }
    if (disc_adj_fem) {
      delete solver_container[val_iInst][iMGlevel][ADJFEA_SOL];
    }
    
    delete [] solver_container[val_iInst][iMGlevel];
  }
  
  delete [] solver_container[val_iInst];

}

void CErrorEstimationDriver::Solver_Deletion(CSolver ****solver_container,
                                    CConfig *config, unsigned short val_iInst) {
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat_fvm, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_fem, disc_adj_heat,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj        = false;
  fem              = false;  disc_adj_fem    = false;
  heat_fvm        = false;   disc_adj_heat   = false;
  transition       = false;
  template_solver  = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : euler = true; break;
    case FEM_NAVIER_STOKES:
    case FEM_LES: ns = true; break;
    case FEM_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent){
    switch (config->GetKind_Turb_Model()) {
    case SA:     spalart_allmaras = true;     break;
    case SA_NEG: neg_spalart_allmaras = true; break;
    case SST:    menter_sst = true;           break;
    case SA_E: e_spalart_allmaras = true; break;
    case SA_COMP: comp_spalart_allmaras = true; break;
    case SA_E_COMP: e_comp_spalart_allmaras = true; break;
    }
  }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
    
    /*--- DeAllocate solution for a template problem ---*/
    
  if (template_solver) {
    delete solver_container[val_iInst][MESH_0][TEMPLATE_SOL];
  }

    /*--- DeAllocate solution for adjoint problem ---*/

  if (adj_euler || adj_ns || disc_adj) {
    delete solver_container[val_iInst][MESH_0][ADJFLOW_SOL];
    if (disc_adj_turb || adj_turb) {
      delete solver_container[val_iInst][MESH_0][ADJTURB_SOL];
    }
    if (heat_fvm) {
      delete solver_container[val_iInst][MESH_0][ADJHEAT_SOL];
    }
  }

  if (disc_adj_heat) {
    delete solver_container[val_iInst][MESH_0][ADJHEAT_SOL];
  }

    /*--- DeAllocate solution for direct problem ---*/

  if (euler || ns) {
    delete solver_container[val_iInst][MESH_0][FLOW_SOL];
  }

  if (turbulent) {
    if (spalart_allmaras || neg_spalart_allmaras || menter_sst || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
      delete solver_container[val_iInst][MESH_0][TURB_SOL];
    }
    if (transition) {
      delete solver_container[val_iInst][MESH_0][TRANS_SOL];
    }
  }
  if (heat_fvm) {
    delete solver_container[val_iInst][MESH_0][HEAT_SOL];
  }
  if (fem) {
    delete solver_container[val_iInst][MESH_0][FEA_SOL];
  }
  if (disc_adj_fem) {
    delete solver_container[val_iInst][MESH_0][ADJFEA_SOL];
  }

  delete [] solver_container[val_iInst][MESH_0];

  delete [] solver_container[val_iInst];

}
