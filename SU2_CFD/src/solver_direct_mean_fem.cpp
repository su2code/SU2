/*!
 * \file solution_direct_mean_fem.cpp
 * \brief Main subroutines for solving finite element flow problems (Euler, Navier-Stokes, etc.).
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 4.0.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/solver_structure.hpp"

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(void) : CSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  
  /*--- Surface-based array initialization ---*/
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL; Surface_CSideForce_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL; Surface_CSideForce = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;
}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  /*--- Determine the restart information. ---*/
  const bool restart = (config->GetRestart() || config->GetRestart_Flow());
  string filename = config->GetSolution_FlowFileName();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL; CEff_Inv = NULL;
  CMx_Inv = NULL;   CMy_Inv = NULL;   CMz_Inv = NULL;
  CFx_Inv = NULL;   CFy_Inv = NULL;   CFz_Inv = NULL;
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL; Surface_CSideForce_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL;   Surface_CFy_Inv = NULL;   Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL;   Surface_CMy_Inv = NULL;   Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL; Surface_CSideForce = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL;   Surface_CFy = NULL;   Surface_CFz = NULL;
  Surface_CMx = NULL;   Surface_CMy = NULL;   Surface_CMz = NULL;
  
  Cauchy_Serie = NULL;
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure. ---*/
  
  nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  const bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);

  if( compressible ) nVar = nDim + 2;
  else               nVar = nDim + 1;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  nMeshPoints = DGGeometry->GetNMeshPoints();
  meshPoints  = DGGeometry->GetMeshPoints();

  nMatchingInternalFaces = DGGeometry->GetNMatchingFaces();
  matchingInternalFaces  = DGGeometry->GetMatchingFaces();

  boundaries = DGGeometry->GetBoundaries();

  standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();
  standardElementsSol      = DGGeometry->GetStandardElementsSol();
  standardMatchingFacesSol = DGGeometry->GetStandardMatchingFacesSol();
  
  /*--- Perform the non-dimensionalization for the flow equations using the
        specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  //LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  //LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Non-dimensional coefficients ---*/
  
  CDrag_Inv         = new su2double[nMarker];
  CLift_Inv         = new su2double[nMarker];
  CSideForce_Inv    = new su2double[nMarker];
  CMx_Inv           = new su2double[nMarker];
  CMy_Inv           = new su2double[nMarker];
  CMz_Inv           = new su2double[nMarker];
  CEff_Inv          = new su2double[nMarker];
  CFx_Inv           = new su2double[nMarker];
  CFy_Inv           = new su2double[nMarker];
  CFz_Inv           = new su2double[nMarker];
  
  Surface_CLift_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CLift          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag          = new su2double[config->GetnMarker_Monitoring()];
  
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag   = 0.0;	Total_CLift        = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy          = 0.0;  Total_CMz          = 0.0;
  Total_CFx     = 0.0;	Total_CFy          = 0.0;  Total_CFz          = 0.0;
  Total_CEff    = 0.0;
  
  /*--- Read farfield conditions ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();

  /*--- Set the conservative variables of the free-stream. ---*/

  vector<su2double> ConsVarFreeStream(nVar);
  if( compressible ) {
    ConsVarFreeStream[0] = Density_Inf;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      ConsVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim];
    ConsVarFreeStream[nVar-1] = Density_Inf*Energy_Inf;
  }
  else {
    ConsVarFreeStream[0] = Pressure_Inf;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      ConsVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim];
  }

  /*--- Determine the total number of DOFs stored on this rank. Allocate the memory
        to store the conservative variables and set the pointers for this solution
        in the elements. ---*/

  nDOFsLocOwned = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) nDOFsLocOwned += volElem[i].nDOFsSol;

  nDOFsLocTot = nDOFsLocOwned;
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) nDOFsLocTot += volElem[i].nDOFsSol;

  VecSolDOFs.resize(nVar*nDOFsLocTot);

  nDOFsLocTot = 0;
  for(unsigned long i=0; i<nVolElemTot; ++i) {
    volElem[i].solDOFs = VecSolDOFs.data() + nVar*nDOFsLocTot;

    nDOFsLocTot += volElem[i].nDOFsSol;
  }
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
  if (!restart || (iMesh != MESH_0)) {
    
    /*--- Start the solution from the free-stream state ---*/

    unsigned long ii = 0;
    for(unsigned long i=0; i<nDOFsLocTot; ++i) {
      for(unsigned short j=0; j<nVar; ++j, ++ii) {
        VecSolDOFs[ii] = ConsVarFreeStream[j];
      }
    }
    
  } else {
    
    /*--- Open the restart file, throw an error if this fails. ---*/

    ifstream restart_file;
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }

    /*--- Create the map from the global DOF ID to the local index. ---*/
    map<unsigned long, unsigned long> mapGlobal2Local;

    unsigned long ii = 0;
    for(unsigned long i=0; i<nVolElemOwned; ++i) {
      for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
        mapGlobal2Local[volElem[i].offsetDOFsSolGlobal+j] = ii;
      }
    }
    
    /*--- The first line is the header ---*/
    string text_line;
    getline (restart_file, text_line);
    
    /*--- Read all lines in the restart file ---*/
    unsigned long iDOF_Global = 0, nDOF_Read = 0;
    while (getline (restart_file, text_line)) {
      istringstream point_line(text_line);

      /*--- Check if this DOF must be stored on this rank. ---*/
      map<unsigned long, unsigned long>::const_iterator MI;
      MI = mapGlobal2Local.find(iDOF_Global);
      if(MI != mapGlobal2Local.end()) {

        /*--- This DOF must be stored on this rank. Retrieve the local index
              and read the data from file. ---*/
        const unsigned long iDOF_Local = nVar*MI->second;

        unsigned long index;
        point_line >> index;
        for(unsigned short i=0; i<nDim; ++i) {
          su2double dull_val;
          point_line >> dull_val;
        }

        for(unsigned short i=0; i<nVar; ++i) {
          point_line >> VecSolDOFs[iDOF_Local+i];
        }

        /*--- Update the local counter nDOF_Read. ---*/
        ++nDOF_Read;
      }

      /*--- Update the counter iDOF_Global. ---*/
      ++iDOF_Global;
    }
    
    /*--- Detect a wrong solution file ---*/

    unsigned short rbuf_NotMatching = 0;
    if(nDOF_Read < nDOFsLocOwned) rbuf_NotMatching = 1;
    
#ifdef HAVE_MPI
    unsigned short sbuf_NotMatching = rbuf_NotMatching;
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  
  unsigned long nBadDOFs = 0;

  if( compressible ) {

    for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
      const unsigned long ii = nVar*i;
      su2double Density = VecSolDOFs[ii];

      su2double Velocity2 = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel = VecSolDOFs[ii+iDim]/Density;
        Velocity2 += vel*vel;
      }

      su2double StaticEnergy = VecSolDOFs[ii+nDim+1]/Density - 0.5*Velocity2;
  
      FluidModel->SetTDState_rhoe(Density, StaticEnergy);
      su2double Pressure = FluidModel->GetPressure();
      su2double Temperature = FluidModel->GetTemperature();
    
      /*--- Use the values at the infinity if the state is not physical. ---*/
    
      if((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
        for(unsigned short j=0; j<nVar; ++j) {
          VecSolDOFs[ii+j] = ConsVarFreeStream[j];
        }

        ++nBadDOFs;
      }
    }
  }
    
  /*--- Warning message about non-physical points ---*/
    
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double nBadDOFsLoc = nBadDOFs;
    SU2_MPI::Reduce(&nBadDOFsLoc, &nBadDOFs, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#endif

    if((rank == MASTER_NODE) && (nBadDOFs != 0))
      cout << "Warning. The original solution contains "<< nBadDOFs << " DOFs that are not physical." << endl;
  }

  /*--- Set up the persistent communication for the conservative variables. ---*/

  Prepare_MPI_Communication(DGGeometry, config);
  
  /*--- Perform the MPI communication of the solution including the possible
        self communication for the periodic data. Correct for rotational
        periodicity afterwards, if needed. ---*/
#ifdef HAVE_MPI
  if( nCommRequests ) {
    MPI_Startall(nCommRequests, commRequests.data());
    SU2_MPI::Waitall(nCommRequests, commRequests.data(), MPI_STATUSES_IGNORE);
  }
#endif

  SelfCommunication();
  CorrectForRotationalPeriodicity();
}

CFEM_DG_EulerSolver::~CFEM_DG_EulerSolver(void) {

  /*--- Array deallocation ---*/
  if (CDrag_Inv != NULL)         delete [] CDrag_Inv;
  if (CLift_Inv != NULL)         delete [] CLift_Inv;
  if (CSideForce_Inv != NULL)    delete [] CSideForce_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (Surface_CLift_Inv != NULL) delete[] Surface_CLift_Inv;
  if (Surface_CDrag_Inv != NULL) delete[] Surface_CDrag_Inv;
  if (Surface_CSideForce_Inv != NULL) delete[] Surface_CSideForce_Inv;
  if (Surface_CEff_Inv != NULL) delete[] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  if (Surface_CLift != NULL)    delete [] Surface_CLift;
  if (Surface_CDrag != NULL)    delete [] Surface_CDrag;
  if (Surface_CSideForce != NULL) delete [] Surface_CSideForce;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;

#ifdef HAVE_MPI

  /*--- Release the memory of the persistent communication and the derived
        data types. ---*/
  for(int i=0; i<nCommRequests; ++i) MPI_Request_free(&commRequests[i]);
  for(int i=0; i<nCommRequests; ++i) MPI_Type_free(&commTypes[i]);

#endif
}

void CFEM_DG_EulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
  
  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
  Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;
  
  unsigned short iDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables ---*/
  
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool standard_air       = (config->GetKind_FluidModel() == STANDARD_AIR);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);
  
  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  
  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();
  
  switch (config->GetKind_FluidModel()) {
      
    case STANDARD_AIR:
      
      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);
      
      FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;
      
    case IDEAL_GAS:
      
      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;
      
    case VW_GAS:
      
      FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                       config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;
      
    case PR_GAS:
      
      FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                     config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;
      
  }
  
  Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();
  
  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  
  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }
  
  /*--- Compute the modulus of the free stream velocity ---*/
  
  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);
  
  /*--- Viscous initialization ---*/
  
  if (viscous) {
    
    /*--- Reynolds based initialization ---*/
    
    if (reynolds_init) {
      
      /*--- First, check if there is mesh motion. If yes, use the Mach
       number relative to the body to initialize the flow. ---*/
      
      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;
      
      /*--- Change of measurement system, hard coded value working only with STANDAR AIR model ---*/
      
      if (standard_air) {
        if (config->GetSystemMeasurements() == SI) {
          config->SetMu_RefND(1.716E-5);
          config->SetMu_SND(110.4);
          config->SetMu_Temperature_RefND(273.15);
        }
        if (config->GetSystemMeasurements() == US) {
          config->SetMu_RefND(3.62E-7);
          config->SetMu_SND(198.72);
          config->SetMu_Temperature_RefND(518.7);
        }
      }
      
      /*--- For viscous flows, pressure will be computed from a density
       that is found from the Reynolds number. The viscosity is computed
       from the dimensional version of Sutherland's law ---*/
      
      FluidModel->SetLaminarViscosityModel(config);
      
      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      
      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = FluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
      
    }
    
    /*--- Thermodynamics quantities based initialization ---*/
    
    else {
      
      FluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
      
    }
    
    /*--- Turbulence kinetic energy ---*/
    
    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
    
  }
  else {
    
    /*--- For inviscid flow, energy is calculated from the specified
     FreeStream quantities using the proper gas law. ---*/
    
    Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
    
  }
  
  /*-- Compute the freestream energy. ---*/
  
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);
  
  /*--- Compute non dimensional quantities. By definition,
   Lref is one because we have converted the grid to meters. ---*/
  
  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);
  
  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);         config->SetFroude(Froude);
  
  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }
  
  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  
  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);
  
  
  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);
  
  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);
  
  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);
  
  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);
  
  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);
  
  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);
  
  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/
  
  switch (config->GetKind_FluidModel()) {
      
    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
  }
  
  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  
  if (viscous) {
    
    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_ConstantND()/Viscosity_Ref);
    
    /*--- Sutherland's model ---*/
    
    config->SetMu_RefND(config->GetMu_RefND()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_SND()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_RefND()/config->GetTemperature_Ref());
    
    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_ConstantND()/Conductivity_Ref);
    
    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);
    
  }
  
  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);
  
  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);
  
  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);
  
  /*--- Write output to the console if this is the master node and first domain ---*/
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {
    
    cout.precision(6);
    
    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }
    
    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;
    
    cout <<"-- Input conditions:"<< endl;
    
    switch (config->GetKind_FluidModel()) {
        
      case STANDARD_AIR:
        cout << "Fluid Model: STANDARD_AIR "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant();
        if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;
        
      case IDEAL_GAS:
        cout << "Fluid Model: IDEAL_GAS "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;
        
      case VW_GAS:
        cout << "Fluid Model: Van der Waals "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;
        
      case PR_GAS:
        cout << "Fluid Model: Peng-Robinson "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;
        
    }
    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {
          
        case CONSTANT_VISCOSITY:
          cout << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
          cout << "Laminar Viscosity: " << config->GetMu_ConstantND()*Viscosity_Ref;
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          break;
          
        case SUTHERLAND:
          cout << "Viscosity Model: SUTHERLAND "<< endl;
          cout << "Ref. Laminar Viscosity: " << config->GetMu_RefND()*Viscosity_Ref;
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Ref. Temperature: " << config->GetMu_Temperature_RefND()*config->GetTemperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Sutherland Constant: "<< config->GetMu_SND()*config->GetTemperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          cout << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< endl;
          cout << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< endl;
          break;
          
      }
      switch (config->GetKind_ConductivityModel()) {
          
        case CONSTANT_PRANDTL:
          cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
          cout << "Prandtl: " << config->GetPrandtl_Lam()<< endl;
          break;
          
        case CONSTANT_CONDUCTIVITY:
          cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
          cout << "Molecular Conductivity: " << config->GetKt_ConstantND()*Conductivity_Ref<< " W/m^2.K." << endl;
          cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
          break;
          
      }
    }
    
    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    
    cout << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";
    
    cout << "Magnitude: "	<< config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;
    
    cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
        cout << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " 1/s." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " 1/s." << endl;
      }
    }
    
    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }
    
    /*--- Print out reference values. ---*/
    
    cout <<"-- Reference values:"<< endl;
    
    cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
    
    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Reference temperature: " << config->GetTemperature_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    
    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;
    
    cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    
    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      cout << "Reference conductivity: " << config->GetConductivity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " W/m^2.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf/ft.s.R." << endl;
    }
    
    
    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;
    
    /*--- Print out resulting non-dim values here. ---*/
    
    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft." << endl;
    }
    
    cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
    cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;
    
    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;
    
    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "	 << config->GetModVel_FreeStreamND() << endl;
    
    cout << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << endl;
        cout << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << endl;
      }
    }
    
    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }
    
    cout << endl;
    
  }
}

void CFEM_DG_EulerSolver::Prepare_MPI_Communication(const CMeshFEM *FEMGeometry,
                                                    CConfig        *config) {

  /*--- Get the communication information from DG_Geometry. Note that for a
        FEM DG discretization the communication entities of FEMGeometry contain
        the volume elements. ---*/
  const vector<int>                    &ranksComm       = FEMGeometry->GetRanksComm();
  const vector<vector<unsigned long> > &elementsSend    = FEMGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsReceive = FEMGeometry->GetEntitiesReceive();

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Find out whether or not self communicationn is present.    ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the rank inside MPI_COMM_WORLD. */
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* Loop over the entries of ranksComm and check for self communication. */
  for(unsigned long i=0; i<ranksComm.size(); ++i) {
    if(ranksComm[i] == rank) {

      /* The send and receive elements for this entry in elementsSend and
         elementsReceive correspond to self communication. Copy the data. */
      elementsReceiveSelfComm = elementsReceive[i];
      elementsSendSelfComm    = elementsSend[i];
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Set up the persistent MPI communication.                   ---*/
  /*--------------------------------------------------------------------------*/

#ifdef HAVE_MPI

  /*--- Determine the number of ranks with which communication takes place,
        multiply this number by two to obtain the number of communication
        requests (send and receive) and allocate the necessary memory. ---*/
  nCommRequests = ranksComm.size();
  if( elementsReceiveSelfComm.size() ) --nCommRequests;
  nCommRequests *= 2;

  commRequests.resize(nCommRequests);
  commTypes.resize(nCommRequests);

  /*--- Loop over the ranks for which communication takes place and exclude
        the self communication. ---*/
  unsigned int nn = 0;
  for(unsigned long i=0; i<ranksComm.size(); ++i) {
    if(ranksComm[i] != rank) {

      /*--- Determine the derived data type for sending the data. ---*/
      const unsigned int nElemSend = elementsSend[i].size();
      vector<int> blockLen(nElemSend), displ(nElemSend);

      for(unsigned int j=0; j<nElemSend; ++j) {
        const unsigned long jj = elementsSend[i][j];
        blockLen[j] = nVar*volElem[jj].nDOFsSol;
        displ[j]    = nVar*volElem[jj].offsetDOFsSolLocal;
      }

      MPI_Type_indexed(nElemSend, blockLen.data(), displ.data(),
                       MPI_DOUBLE, &commTypes[nn]);
      MPI_Type_commit(&commTypes[nn]);

      /* Create the communication request for this send operation. Afterwards
         update the counter nn for the next request. */
      MPI_Send_init(VecSolDOFs.data(), 1, commTypes[nn], ranksComm[i],
                    ranksComm[i], MPI_COMM_WORLD, &commRequests[nn]);
      ++nn;

      /*--- Determine the derived data type for receiving the data. ---*/
      const unsigned int nElemRecv = elementsReceive[i].size();
      blockLen.resize(nElemRecv);
      displ.resize(nElemRecv);

      for(unsigned int j=0; j<nElemRecv; ++j) {
        const unsigned long jj = elementsReceive[i][j];
        blockLen[j] = nVar*volElem[jj].nDOFsSol;
        displ[j]    = nVar*volElem[jj].offsetDOFsSolLocal;
      }

      MPI_Type_indexed(nElemRecv, blockLen.data(), displ.data(),
                       MPI_DOUBLE, &commTypes[nn]);
      MPI_Type_commit(&commTypes[nn]);

      /* Create the communication request for this receive operation. Afterwards
         update the counter nn for the next request. */
      MPI_Recv_init(VecSolDOFs.data(), 1, commTypes[nn], ranksComm[i],
                    rank, MPI_COMM_WORLD, &commRequests[nn]);
      ++nn;
    }
  }

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Store the information for the rotational periodic          ---*/
  /*---         corrections for the halo elements.                         ---*/
  /*--------------------------------------------------------------------------*/

  /* Get the data for the rotational periodic halos from FEMGeometry. */
  vector<unsigned short> markersRotationalPeriodicity = FEMGeometry->GetRotPerMarkers();

  halosRotationalPeriodicity = FEMGeometry->GetRotPerHalos();

  /* Determine the number of rotational periodic transformations and allocate
     the memory to store the corresponding rotation matrices. */
  const unsigned short nRotPerMarkers = markersRotationalPeriodicity.size();
  rotationMatricesPeriodicity.resize(9*nRotPerMarkers);

  /* Loop over the rotational periodic transformations and store the rotation
     matrices in rotationMatricesPeriodicity. */
  unsigned int ii = 0;
  for(unsigned short i=0; i<nRotPerMarkers; ++i) {

   /* Get the rotation angles from config for this marker. */
   const unsigned short pInd = markersRotationalPeriodicity[i];
   su2double *angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(pInd));

    /*--- Determine the rotation matrix from the donor to the halo elements.
          This is the transpose of the rotation matrix from the halo to the
          donor elements, which is stored in periodic angles of the marker. ---*/
    const su2double theta = angles[0], phi = angles[1], psi = angles[2];

    const su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
    const su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

    rotationMatricesPeriodicity[ii++] =  cosPhi*cosPsi;
    rotationMatricesPeriodicity[ii++] =  cosPhi*sinPsi;
    rotationMatricesPeriodicity[ii++] = -sinPhi;

    rotationMatricesPeriodicity[ii++] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotationMatricesPeriodicity[ii++] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotationMatricesPeriodicity[ii++] = sinTheta*cosPhi;

    rotationMatricesPeriodicity[ii++] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotationMatricesPeriodicity[ii++] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotationMatricesPeriodicity[ii++] = cosTheta*cosPhi;
  }
}

void CFEM_DG_EulerSolver::CorrectForRotationalPeriodicity(void) {

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii;
  const unsigned short nRotPerMarkers = halosRotationalPeriodicity.size();
  for(unsigned short k=0; k<nRotPerMarkers; ++k) {

    /* Easier storage of the rotational matrix. */
    su2double rotMatrix[3][3];

    rotMatrix[0][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[0][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[0][2] = rotationMatricesPeriodicity[ii++];

    rotMatrix[1][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][2] = rotationMatricesPeriodicity[ii++];

    rotMatrix[2][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][2] = rotationMatricesPeriodicity[ii++];

    /* Determine the number of elements for this transformation and
       loop over them. */
    const unsigned long nHaloElem = halosRotationalPeriodicity[k].size();
    for(unsigned long j=0; j<nHaloElem; ++j) {

      /* Easier storage of the halo index and loop over its DOFs. */
      const unsigned long ind = halosRotationalPeriodicity[k][j];
      for(unsigned short i=0; i<volElem[ind].nDOFsSol; ++i) {

        /* Determine the pointer in solDOFs where the solution of this DOF
           is stored and copy the momentum variables. Note that a rotational
           correction can only take place for a 3D simulation. */
        su2double *sol = volElem[ind].solDOFs + i*nVar;

        const su2double ru = sol[1], rv = sol[2], rw = sol[3];

        /* Correct the momentum variables. */
        sol[1] = rotMatrix[0][0]*ru + rotMatrix[0][1]*rv + rotMatrix[0][2]*rw;
        sol[2] = rotMatrix[1][0]*ru + rotMatrix[1][1]*rv + rotMatrix[1][2]*rw;
        sol[3] = rotMatrix[2][0]*ru + rotMatrix[2][1]*rv + rotMatrix[2][2]*rw;
      }
    }
  }
}

void CFEM_DG_EulerSolver::SelfCommunication(void) {

  /* Easier storage of the number of variables times the size of su2double. */
  const unsigned long sizeEntity = nVar*sizeof(su2double);

  /* Loop over the number of elements involved and copy the data of the DOFs. */
  const unsigned long nSelfElements = elementsSendSelfComm.size();
  for(unsigned long i=0; i<nSelfElements; ++i) {
    const unsigned long indS   = elementsSendSelfComm[i];
    const unsigned long indR   = elementsReceiveSelfComm[i];
    const unsigned long nBytes = volElem[indS].nDOFsSol * sizeEntity;

    memcpy(volElem[indR].solDOFs, volElem[indS].solDOFs, nBytes);
  }
}

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply any problem-dependent initial conditions here. ---*/
  
  if (ExtIter == 0 && rank == MASTER_NODE) cout << " Setting initial condition." << endl;
  
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long ErrorCounter = 0;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Perform any preprocessing here, e.g., compute gradients, limiters, special dissipation, etc. ---*/
  
  if (rank == MASTER_NODE) cout << " Preprocessing." << endl;
  
  /*--- Collect the number of non-physical points for this iteration. ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CFEM_DG_EulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) { }

unsigned long CFEM_DG_EulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    RightSol = node[iPoint]->SetPrimVar_Compressible(FluidModel);
    node[iPoint]->SetSecondaryVar_Compressible(FluidModel);
    
    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Compute and store time step here. ---*/
  
  if (rank == MASTER_NODE) cout << " Setting time step." << endl;
  
}

void CFEM_DG_EulerSolver::Convective_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                              CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Compute finite element convective residual and store here. ---*/
  
  if (rank == MASTER_NODE) cout << " Computing convective residual." << endl;
  
  // numerics holds the CDiscGalerkin class for compute residual here
  
  if ((config->GetRiemann_Solver_FEM() == ROE) && (rank == MASTER_NODE))
      cout << " Using a Roe Riemann solver for the DG method." << endl;
  
  su2double quad_fact_straight = config->GetQuadrature_Factor_Straight();
  if (rank == MASTER_NODE) cout << " Quad factor straight: " << quad_fact_straight << endl;
  
  su2double quad_fact_curved = config->GetQuadrature_Factor_Curved();
  if (rank == MASTER_NODE) cout << " Quad factor curved: " << quad_fact_curved << endl;
  
}

void CFEM_DG_EulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                       CConfig *config, unsigned short iMesh) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Compute source residual and store here. ---*/
  
  if (rank == MASTER_NODE) cout << " Computing source residual." << endl;
  
  // numerics holds the CSourceFEM class for source resid here
  
}

void CFEM_DG_EulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) cout << " Computing inviscid forces." << endl;
  
}

void CFEM_DG_EulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                            CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CFEM_DG_EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply flow tangency boundary condition here. ---*/
  
  if (rank == MASTER_NODE) cout << " Euler wall BC." << endl;
  
}

void CFEM_DG_EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply farfield boundary condition here. ---*/
  
  if (rank == MASTER_NODE) cout << " Far-field BC." << endl;
  
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler residual ---*/
  
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
  
}


CFEM_DG_NSSolver::CFEM_DG_NSSolver(void) : CFEM_DG_EulerSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Visc = NULL; CLift_Visc = NULL; CSideForce_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;
  
  /*--- Surface-based array initialization ---*/
  
  Surface_CLift_Visc = NULL; Surface_CDrag_Visc = NULL; Surface_CSideForce_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  
}

CFEM_DG_NSSolver::CFEM_DG_NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
 : CFEM_DG_EulerSolver(geometry, config, iMesh) {
  
  /*--- Array initialization ---*/
  
  CDrag_Visc = NULL; CLift_Visc = NULL; CSideForce_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  Surface_CLift_Visc = NULL; Surface_CDrag_Visc = NULL; Surface_CSideForce_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  
  ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction = NULL;    Cauchy_Serie = NULL;

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  //LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  //LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Non dimensional coefficients ---*/

  ForceViscous     = new su2double[3];
  MomentViscous    = new su2double[3];
  CDrag_Visc       = new su2double[nMarker];
  CLift_Visc       = new su2double[nMarker];
  CSideForce_Visc  = new su2double[nMarker];
  CMx_Visc         = new su2double[nMarker];
  CMy_Visc         = new su2double[nMarker];
  CMz_Visc         = new su2double[nMarker];
  CEff_Visc        = new su2double[nMarker];
  CFx_Visc         = new su2double[nMarker];
  CFy_Visc         = new su2double[nMarker];
  CFz_Visc         = new su2double[nMarker];
  
  Surface_CLift_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag   = 0.0;	Total_CLift        = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy          = 0.0;  Total_CMz          = 0.0;
  Total_CEff    = 0.0;
  Total_CFx     = 0.0;	Total_CFy          = 0.0;  Total_CFz          = 0.0;
  
  /*--- Read farfield conditions from config ---*/
  
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();
}

CFEM_DG_NSSolver::~CFEM_DG_NSSolver(void) {
  unsigned short iMarker;
  
  if (CDrag_Visc != NULL)       delete [] CDrag_Visc;
  if (CLift_Visc != NULL)       delete [] CLift_Visc;
  if (CSideForce_Visc != NULL)  delete [] CSideForce_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;
  
  
  if (Surface_CLift_Visc != NULL)      delete [] Surface_CLift_Visc;
  if (Surface_CDrag_Visc != NULL)      delete [] Surface_CDrag_Visc;
  if (Surface_CSideForce_Visc != NULL) delete [] Surface_CSideForce_Visc;
  if (Surface_CEff_Visc != NULL)       delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)        delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)        delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)        delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)        delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)        delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)        delete [] Surface_CMz_Visc;
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;
  
  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }
  
}

void CFEM_DG_NSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  bool limiter_visc = config->GetViscous_Limiter_Flow();
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  StrainMag_Max = 0.0, Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity(limiter_visc);
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag(limiter_visc);
    
    StrainMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }
  
  /*--- Collect the number of non-physical points for this iteration. ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    su2double MyOmega_Max = Omega_Max; Omega_Max = 0.0;
    su2double MyStrainMag_Max = StrainMag_Max; StrainMag_Max = 0.0;
    
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    if (iMesh == MESH_0) {
      config->SetNonphysical_Points(ErrorCounter);
      solver_container[FLOW_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[FLOW_SOL]->SetOmega_Max(Omega_Max);
    }
    
  }
  
}

void CFEM_DG_NSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) cout << " Computing viscous forces." << endl;
}

void CFEM_DG_NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  
}

void CFEM_DG_NSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  // numerics holds the CFEMVisc class for viscous resid here
  
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply adiabatic / heat flux wall boundary condition here. ---*/
  
  if (rank == MASTER_NODE) cout << " Heat flux wall BC." << endl;
  
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply isothermal wall boundary condition here. ---*/
  
  if (rank == MASTER_NODE) cout << " Isothermal wall BC." << endl;
  
}
