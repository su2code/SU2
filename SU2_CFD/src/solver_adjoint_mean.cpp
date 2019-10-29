/*!
 * \file solution_adjoint_mean.cpp
 * \brief Main subroutines for solving adjoint problems (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon, H. Kline
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

#include "../include/solver_structure.hpp"
#include "../include/variables/CAdjEulerVariable.hpp"
#include "../include/variables/CAdjNSVariable.hpp"

CAdjEulerSolver::CAdjEulerSolver(void) : CSolver() {
  
  /*--- Array initialization ---*/
  Phi_Inf = NULL;
  Sens_Mach = NULL;
  Sens_AoA = NULL;
  Sens_Geo = NULL;
  Sens_Press = NULL;
  Sens_Temp = NULL;
  Sens_BPress = NULL;
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  Jacobian_Axisymmetric = NULL;
  CSensitivity = NULL;
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  DonorAdjVar = NULL;
  DonorGlobalIndex = NULL;

}

CAdjEulerSolver::CAdjEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  unsigned long iPoint, iVertex;
  string text_line, mesh_filename;
  unsigned short iDim, iVar, iMarker, nLineLets;
  ifstream restart_file;
  string filename, AdjExt;
  su2double myArea_Monitored, Area, *Normal;

  adjoint = true;

  bool restart  = config->GetRestart();

  bool axisymmetric = config->GetAxisymmetric();
  
  su2double RefArea    = config->GetRefArea();
  su2double RefDensity  = config->GetDensity_FreeStreamND();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double Mach_Motion     = config->GetMach_Motion();
  su2double RefVel2, Mach2Vel, Weight_ObjFunc, factor;
  su2double *Velocity_Inf;
  string Marker_Tag, Monitoring_Tag;
  unsigned short iMarker_Monitoring, jMarker, ObjFunc;
  bool grid_movement  = config->GetGrid_Movement();

  /*--- Array initialization ---*/
  
  Phi_Inf = NULL;
  Sens_Mach = NULL;
  Sens_AoA = NULL;
  Sens_Geo = NULL;
  Sens_Press = NULL;
  Sens_Temp = NULL;
  Sens_BPress = NULL;
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  Jacobian_Axisymmetric = NULL;
  CSensitivity = NULL;
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  DonorAdjVar = NULL;
  DonorGlobalIndex = NULL;

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constans in the solver structure ---*/
  nDim = geometry->GetnDim();
  nMarker = config->GetnMarker_All();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  nVar = nDim + 2;
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  Residual_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
  Res_Conv_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]   = 0.0;
  Solution_j = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary arrays related to the flow solution ---*/
  FlowPrimVar_i = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_i[iVar] = 0.0;
  FlowPrimVar_j = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
    Jacobian_ii = new su2double* [nVar];
    Jacobian_ij = new su2double* [nVar];
    Jacobian_ji = new su2double* [nVar];
    Jacobian_jj = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_ii[iVar] = new su2double [nVar];
      Jacobian_ij[iVar] = new su2double [nVar];
      Jacobian_ji[iVar] = new su2double [nVar];
      Jacobian_jj[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    if (axisymmetric) {
      Jacobian_Axisymmetric = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_Axisymmetric[iVar] = new su2double [nVar];
    }
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Sensitivity definition and coefficient in all the markers ---*/
  CSensitivity = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]];
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  DonorAdjVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorAdjVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorAdjVar[iMarker][iVertex] = new su2double [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        DonorAdjVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }
  
  Sens_Geo  = new su2double[nMarker];
  Sens_Mach = new su2double[nMarker];
  Sens_AoA  = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];
  Sens_BPress = new su2double[nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker]  = 0.0;
    Sens_Mach[iMarker] = 0.0;
    Sens_AoA[iMarker]  = 0.0;
    Sens_Press[iMarker] = 0.0;
    Sens_Temp[iMarker]  = 0.0;
    Sens_BPress[iMarker] = 0.0;
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
      CSensitivity[iMarker][iVertex] = 0.0;
  }

  /*--- Adjoint flow at the inifinity, initialization stuff ---*/
  PsiRho_Inf = 0.0; PsiE_Inf   = 0.0;
  Phi_Inf    = new su2double [nDim];
  Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
  if (nDim == 3) Phi_Inf[2] = 0.0;
  
  /*--- If outflow objective, nonzero initialization ---*/
  if ((config->GetKind_ObjFunc() == SURFACE_TOTAL_PRESSURE)) {
    su2double SoundSpeed,*vel_inf,R,vel2,vel;
    R = config->GetGas_ConstantND();
    vel_inf = config->GetVelocity_FreeStreamND();
    vel2=0;
    for (iDim=0; iDim<nDim; iDim++)
      vel2 +=vel_inf[iDim]*vel_inf[iDim];
    vel = pow(vel2,0.5);
    SoundSpeed= pow(Gamma*config->GetTemperature_FreeStreamND()*R, 0.5);
    PsiE_Inf = Gamma_Minus_One*vel2/(vel2-pow(SoundSpeed,2.0))*0.5/vel;
    PsiRho_Inf += PsiE_Inf*(2*SoundSpeed*SoundSpeed+vel2*Gamma_Minus_One)/(2.0*Gamma_Minus_One);
    // Assumes +x flow direction
    // Assume v.n = |v|, n = -v/|v|

    for (iDim=0; iDim<nDim; iDim++) {
      Phi_Inf[iDim] +=PsiE_Inf*(SoundSpeed*SoundSpeed/Gamma_Minus_One/vel2-1)*vel_inf[iDim];
      // Assumes n in direction of v
      Phi_Inf[iDim]+=vel_inf[iDim]/vel*(0.5);
    }

  }

  /*--- Initialize the solution with the far-field state everywhere. ---*/

  nodes = new CAdjEulerVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Read the restart metadata. ---*/

  if (restart && (iMesh == MESH_0)) {
    mesh_filename = config->GetSolution_AdjFileName();
    filename      = config->GetObjFunc_Extension(mesh_filename);
//    Read_SU2_Restart_Metadata(geometry, config, true, filename);
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;


  /*--- Calculate area monitored for area-averaged-outflow-quantity-based objectives ---*/

  myArea_Monitored = 0.0;
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    if (config->GetKind_ObjFunc(iMarker_Monitoring)==SURFACE_TOTAL_PRESSURE ||
        config->GetKind_ObjFunc(iMarker_Monitoring)==SURFACE_STATIC_PRESSURE) {

      Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
      /*-- Find the marker index ---*/
      iMarker = 0;
      for (jMarker= 0; jMarker < config->GetnMarker_All(); jMarker++) {
        Marker_Tag = config->GetMarker_All_TagBound(jMarker);
        if (Marker_Tag == Monitoring_Tag) {
          iMarker = jMarker;
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            if (geometry->node[iPoint]->GetDomain()) {
              Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
              Area = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                Area += Normal[iDim]*Normal[iDim];
              myArea_Monitored += sqrt (Area);
            }
          }
          break;
        }
      }
    }
  }

#ifdef HAVE_MPI
  Area_Monitored = 0.0;
  SU2_MPI::Allreduce(&myArea_Monitored, &Area_Monitored, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  Area_Monitored = myArea_Monitored;
#endif

  if (config->GetnObj() > 1 && iMesh==MESH_0) {
    if (grid_movement) {
      Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
      Velocity_Inf = config->GetVelocity_FreeStreamND();
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }

    /*--- Objective scaling: a factor must be applied to certain objectives ---*/
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);

        factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);

        ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);
        if ((ObjFunc == INVERSE_DESIGN_HEATFLUX) ||
            (ObjFunc == TOTAL_HEATFLUX) || (ObjFunc == MAXIMUM_HEATFLUX) ||
            (ObjFunc == SURFACE_MASSFLOW) ) factor = 1.0;

       if ((ObjFunc == SURFACE_TOTAL_PRESSURE) || (ObjFunc == SURFACE_STATIC_PRESSURE))
         factor = 1.0/Area_Monitored;

       Weight_ObjFunc = Weight_ObjFunc*factor;
       config->SetWeight_ObjFunc(iMarker_Monitoring, Weight_ObjFunc);
    }
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
}

CAdjEulerSolver::~CAdjEulerSolver(void) {
  unsigned short iVar, iMarker;
  
  if (Phi_Inf != NULL) delete [] Phi_Inf;
  if (Sens_Mach != NULL) delete [] Sens_Mach;
  if (Sens_AoA != NULL) delete [] Sens_AoA;
  if (Sens_Geo != NULL) delete [] Sens_Geo;
  if (Sens_Press != NULL) delete [] Sens_Press;
  if (Sens_Temp != NULL) delete [] Sens_Temp;
  if (Sens_BPress != NULL) delete [] Sens_BPress;
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  
  if (Jacobian_Axisymmetric != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_Axisymmetric[iVar];
    delete [] Jacobian_Axisymmetric;
  }
  
  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CSensitivity[iMarker];
    delete [] CSensitivity;
  }
  
  if (nodes != nullptr) delete nodes;
}

void CAdjEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh, unsigned long Iteration) {

  /*--- Use the flow solution to update the time step
   *    The time step depends on the characteristic velocity, which is the same
   *    for the adjoint and flow solutions, albeit in the opposite direction. ---*/
  solver_container[FLOW_SOL]->SetTime_Step(geometry, solver_container, config, iMesh, Iteration);
}

void CAdjEulerSolver::Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;

#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  

#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_AdjVar          = NULL;
  su2double        *iAdjVar          = new su2double [nVar];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_AdjVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_AdjVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_AdjVar          += nPointTotal_s[iDomain]*(nVar+3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_AdjVar          = new su2double[Buffer_Size_AdjVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = nodes->GetSolution(iPoint,iVar);
            }
            Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+(nVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+(nVar+1)] = su2double(jVertex);
            Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+(nVar+2)]  = su2double(jMarker);
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_AdjVar[PointTotal_Counter*(nVar+3)],
                     nPointTotal_s[iDomain]*(nVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_AdjVar            = new su2double[nPointTotal_s[iDomain]*(nVar+3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nVar+3); iter++)
        Buffer_Receive_AdjVar[iter] = Buffer_Send_AdjVar[PointTotal_Counter*(nVar+3)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+2)]);
        for (iVar = 0; iVar < nVar; iVar++)
          iAdjVar[iVar] = Buffer_Receive_AdjVar[iPoint*(nVar+3)+iVar];
        
        for (iVar = 0; iVar < nVar; iVar++)
          SetDonorAdjVar(iMarker, iVertex, iVar, iAdjVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_AdjVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_AdjVar            = new su2double [nPointTotal_r[iDomain]*nVar];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_AdjVar, nPointTotal_r[iDomain]*(nVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+2)]);
        for (iVar = 0; iVar < nVar; iVar++)
          iAdjVar[iVar] = Buffer_Receive_AdjVar[iPoint*(nVar+3)+iVar];
        
        for (iVar = 0; iVar < nVar; iVar++)
          SetDonorAdjVar(iMarker, iVertex, iVar, iAdjVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_AdjVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_AdjVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iAdjVar;
  
}

void CAdjEulerSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;

#ifdef HAVE_MPI

  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  

#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_AdjVar          = NULL;
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  su2double        *iAdjVar          = new su2double [nVar];
  
  unsigned long Buffer_Size_AdjVar          = 0;
  
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  
  su2double        *Buffer_Receive_AdjVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_AdjVar          += nPointTotal_s[iDomain]*(nVar+3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_AdjVar          = new su2double[Buffer_Size_AdjVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = nodes->GetSolution(iPoint,iVar);
            }
            Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+(nVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+(nVar+1)] = su2double(jVertex);
            Buffer_Send_AdjVar[(nVar+3)*(PointTotal_Counter+iPointTotal)+(nVar+2)]  = su2double(jMarker);
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_AdjVar[PointTotal_Counter*(nVar+3)],
                     nPointTotal_s[iDomain]*(nVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_AdjVar            = new su2double[nPointTotal_s[iDomain]*(nVar+3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nVar+3); iter++)
        Buffer_Receive_AdjVar[iter] = Buffer_Send_AdjVar[PointTotal_Counter*(nVar+3)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+2)]);
        for (iVar = 0; iVar < nVar; iVar++)
          iAdjVar[iVar] = Buffer_Receive_AdjVar[iPoint*(nVar+3)+iVar];
        
        for (iVar = 0; iVar < nVar; iVar++)
          SetDonorAdjVar(iMarker, iVertex, iVar, iAdjVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_AdjVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_AdjVar            = new su2double [nPointTotal_r[iDomain]*(nVar+3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_AdjVar, nPointTotal_r[iDomain]*(nVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_AdjVar[iPoint*(nVar+3)+(nVar+2)]);
        for (iVar = 0; iVar < nVar; iVar++)
          iAdjVar[iVar] = Buffer_Receive_AdjVar[iPoint*(nVar+3)+iVar];
        
        for (iVar = 0; iVar < nVar; iVar++)
          SetDonorAdjVar(iMarker, iVertex, iVar,  iAdjVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_AdjVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_AdjVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iAdjVar;
  
}

void CAdjEulerSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  su2double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, CD, CL, Cp, CpTarget,
  CT, CQ, x_origin, y_origin, z_origin, WDrag, Area, invCD, CLCD2, invCQ, CTRCQ2;
  unsigned short iMarker,jMarker,iMarker_Monitoring, iDim;
  unsigned long iVertex, iPoint;
  string Marker_Tag, Monitoring_Tag;
  su2double Weight_ObjFunc=1.0;
  su2double *ForceProj_Vector2;
  
  su2double Alpha            = (config->GetAoA()*PI_NUMBER)/180.0;
  su2double Beta             = (config->GetAoS()*PI_NUMBER)/180.0;
  su2double RefLength  = config->GetRefLength();
  su2double *RefOriginMoment = config->GetRefOriginMoment(0);
  su2double dCD_dCL          = config->GetdCD_dCL();
  su2double dCMx_dCL         = config->GetdCMx_dCL();
  su2double dCMy_dCL         = config->GetdCMy_dCL();
  su2double dCMz_dCL         = config->GetdCMz_dCL();
  su2double dCD_dCMy         = config->GetdCD_dCMy();
  bool Fixed_CL              = config->GetFixed_CL_Mode();
  bool Fixed_CM              = config->GetFixed_CM_Mode();

  ForceProj_Vector = new su2double[nDim];
  
  /*--- Compute coefficients needed for objective function evaluation. ---*/
  
  CD = solver_container[FLOW_SOL]->GetTotal_CD();
  CL = solver_container[FLOW_SOL]->GetTotal_CL();
  CT = solver_container[FLOW_SOL]->GetTotal_CT();
  CQ = solver_container[FLOW_SOL]->GetTotal_CQ();
  invCD  = 1.0/CD; CLCD2  = CL/(CD*CD);
  invCQ  = 1.0/CQ; CTRCQ2 = CT/(RefLength*CQ*CQ);
  
  x_origin = RefOriginMoment[0]; y_origin = RefOriginMoment[1]; z_origin = RefOriginMoment[2];
  
  /*--- Evaluate the boundary condition coefficients,
   Since there may be more than one objective per marker, first we have to set all Force projection vectors to 0 ---*/
  
  for (iMarker = 0; iMarker<nMarker; iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
           (config->GetMarker_All_Monitoring(iMarker) == YES))
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        for (iDim=0; iDim<nDim; iDim++)
          ForceProj_Vector[iDim]=0.0;
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        nodes->SetForceProj_Vector(iPoint,ForceProj_Vector);
      }
  }

  /*--- Find the matching iMarker ---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
    for (jMarker=0; jMarker<nMarker; jMarker++) {
      Marker_Tag = config->GetMarker_All_TagBound(jMarker);
      if (Monitoring_Tag==Marker_Tag)
        iMarker = jMarker;
    }

    
    if ((iMarker<nMarker) && (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
        (config->GetMarker_All_Monitoring(iMarker) == YES)) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        x = geometry->node[iPoint]->GetCoord(0);
        y = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) z = geometry->node[iPoint]->GetCoord(2);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        ForceProj_Vector2 = nodes->GetForceProj_Vector(iPoint);
        for (iDim=0; iDim<nDim;iDim++)
          ForceProj_Vector[iDim]=ForceProj_Vector2[iDim];

        switch (config->GetKind_ObjFunc(iMarker_Monitoring)) {
          case DRAG_COEFFICIENT :
            if (nDim == 2) {
              
              ForceProj_Vector[0] += Weight_ObjFunc*cos(Alpha);
              ForceProj_Vector[1] += Weight_ObjFunc*sin(Alpha);
              
              /*--- Modification to run at a fixed CL and CM value ---*/
              
              if (Fixed_CL) { ForceProj_Vector[0] += dCD_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= dCD_dCL*Weight_ObjFunc*cos(Alpha); }
              if (Fixed_CM) { ForceProj_Vector[0] -= dCD_dCMy*Weight_ObjFunc*(y - y_origin)/RefLength; ForceProj_Vector[1] += dCD_dCMy*Weight_ObjFunc*(x - x_origin)/RefLength; }

            }
            if (nDim == 3) {
              
              ForceProj_Vector[0] += Weight_ObjFunc*cos(Alpha)*cos(Beta);
              ForceProj_Vector[1] += Weight_ObjFunc*sin(Beta);
              ForceProj_Vector[2] += Weight_ObjFunc*sin(Alpha)*cos(Beta);
              
              /*--- Modification to run at a fixed CL value ---*/
              
              if (Fixed_CL) { ForceProj_Vector[0] += dCD_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= 0.0; ForceProj_Vector[2] -= dCD_dCL*Weight_ObjFunc*cos(Alpha); }
              if (Fixed_CM) { ForceProj_Vector[0] += dCD_dCMy*Weight_ObjFunc*(z - z_origin)/RefLength; ForceProj_Vector[1] += 0.0; ForceProj_Vector[2] -= dCD_dCMy*Weight_ObjFunc*(x - x_origin)/RefLength; }
              
            }
            break;
          case LIFT_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] += -Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] += Weight_ObjFunc*cos(Alpha); }
            if (nDim == 3) { ForceProj_Vector[0] += -Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] += Weight_ObjFunc*cos(Alpha); }
            break;
          case SIDEFORCE_COEFFICIENT :
            if (nDim == 2) { SU2_MPI::Error("This functional is not possible in 2D!!", CURRENT_FUNCTION);}
            if (nDim == 3) { ForceProj_Vector[0] += -Weight_ObjFunc*sin(Beta) * cos(Alpha); ForceProj_Vector[1] += Weight_ObjFunc*cos(Beta); ForceProj_Vector[2] += -Weight_ObjFunc*sin(Beta) * sin(Alpha); }
            break;
          case INVERSE_DESIGN_PRESSURE :
            Cp = solver_container[FLOW_SOL]->GetCPressure(iMarker, iVertex);
            CpTarget = solver_container[FLOW_SOL]->GetCPressureTarget(iMarker, iVertex);
            Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1]);
            if (nDim == 3) Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
            ForceProj_Vector[0] += -Weight_ObjFunc*2.0*(Cp-CpTarget)*Normal[0]/Area; ForceProj_Vector[1] += -Weight_ObjFunc*2.0*(Cp-CpTarget)*Normal[1]/Area;
            if (nDim == 3) ForceProj_Vector[2] += -Weight_ObjFunc*2.0*(Cp-CpTarget)*Normal[2]/Area;
            break;
          case MOMENT_X_COEFFICIENT :
            if (nDim == 2) { SU2_MPI::Error("This functional is not possible in 2D!!", CURRENT_FUNCTION);}
            if (nDim == 3) {
              
              ForceProj_Vector[0] += 0.0;
              ForceProj_Vector[1] += -Weight_ObjFunc*(z - z_origin)/RefLength;
              ForceProj_Vector[2] += Weight_ObjFunc*(y - y_origin)/RefLength;
            
              /*--- Modification to run at a fixed CL value ---*/
              
              if (Fixed_CL) { ForceProj_Vector[0] += dCMx_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= 0.0; ForceProj_Vector[2] -= dCMx_dCL*Weight_ObjFunc*cos(Alpha); }
            
            }
            break;
          case MOMENT_Y_COEFFICIENT :
            if (nDim == 2) { SU2_MPI::Error("This functional is not possible in 2D!!", CURRENT_FUNCTION);}
            if (nDim == 3) {
              
              ForceProj_Vector[0] += Weight_ObjFunc*(z - z_origin)/RefLength;
              ForceProj_Vector[1] += 0.0;
              ForceProj_Vector[2] += -Weight_ObjFunc*(x - x_origin)/RefLength;
            
              /*--- Modification to run at a fixed CL value ---*/
              
              if (Fixed_CL) { ForceProj_Vector[0] += dCMy_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= 0.0; ForceProj_Vector[2] -= dCMy_dCL*Weight_ObjFunc*cos(Alpha); }

            }
            break;
          case MOMENT_Z_COEFFICIENT :
            if (nDim == 2) {
              
              ForceProj_Vector[0] += -Weight_ObjFunc*(y - y_origin)/RefLength;
              ForceProj_Vector[1] += Weight_ObjFunc*(x - x_origin)/RefLength;
            
              /*--- Modification to run at a fixed CL and CM value ---*/
              
              if (Fixed_CL) { ForceProj_Vector[0] += dCMz_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= dCMz_dCL*Weight_ObjFunc*cos(Alpha); }

            }
            if (nDim == 3) {
              
              ForceProj_Vector[0] += -Weight_ObjFunc*(y - y_origin)/RefLength;
              ForceProj_Vector[1] += Weight_ObjFunc*(x - x_origin)/RefLength;
              ForceProj_Vector[2] += 0;
            
              /*--- Modification to run at a fixed CL value ---*/
              
              if (Fixed_CL) { ForceProj_Vector[0] += dCMz_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= 0.0; ForceProj_Vector[2] -= dCMz_dCL*Weight_ObjFunc*cos(Alpha); }

            }
            break;
          case EFFICIENCY :
            if (nDim == 2) { ForceProj_Vector[0] += -Weight_ObjFunc*(invCD*sin(Alpha)+CLCD2*cos(Alpha)); ForceProj_Vector[1] += Weight_ObjFunc*(invCD*cos(Alpha)-CLCD2*sin(Alpha)); }
            if (nDim == 3) { ForceProj_Vector[0] += -Weight_ObjFunc*(invCD*sin(Alpha)+CLCD2*cos(Alpha)*cos(Beta)); ForceProj_Vector[1] += -Weight_ObjFunc*CLCD2*sin(Beta); ForceProj_Vector[2] += Weight_ObjFunc*(invCD*cos(Alpha)-CLCD2*sin(Alpha)*cos(Beta)); }
            break;
          case EQUIVALENT_AREA :
            WDrag = Weight_ObjFunc*config->GetWeightCd();
            if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha)*WDrag; ForceProj_Vector[1] = sin(Alpha)*WDrag; }
            if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta)*WDrag; ForceProj_Vector[1] = sin(Beta)*WDrag; ForceProj_Vector[2] = sin(Alpha)*cos(Beta)*WDrag; }
            break;
          case NEARFIELD_PRESSURE :
            WDrag = Weight_ObjFunc*config->GetWeightCd();
            if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha)*WDrag; ForceProj_Vector[1] = sin(Alpha)*WDrag; }
            if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta)*WDrag; ForceProj_Vector[1] = sin(Beta)*WDrag; ForceProj_Vector[2] = sin(Alpha)*cos(Beta)*WDrag; }
            break;
          case FORCE_X_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] += Weight_ObjFunc*1.0; ForceProj_Vector[1] += 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] += Weight_ObjFunc*1.0; ForceProj_Vector[1] += 0.0; ForceProj_Vector[2] += 0.0; }
            break;
          case FORCE_Y_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] += 0.0; ForceProj_Vector[1] += Weight_ObjFunc*1.0; }
            if (nDim == 3) { ForceProj_Vector[0] += 0.0; ForceProj_Vector[1] += Weight_ObjFunc*1.0; ForceProj_Vector[2] += 0.0; }
            break;
          case FORCE_Z_COEFFICIENT :
            if (nDim == 2) { SU2_MPI::Error("This functional is not possible in 2D!!", CURRENT_FUNCTION);}
            if (nDim == 3) { ForceProj_Vector[0] += 0.0; ForceProj_Vector[1] += 0.0; ForceProj_Vector[2] += Weight_ObjFunc*1.0; }
            break;
          case THRUST_COEFFICIENT :
            if (nDim == 2) { SU2_MPI::Error("This functional is not possible in 2D!!", CURRENT_FUNCTION);}
            if (nDim == 3) { ForceProj_Vector[0] += 0.0; ForceProj_Vector[1] += 0.0; ForceProj_Vector[2] += Weight_ObjFunc*1.0; }
            break;
          case TORQUE_COEFFICIENT :
            if (nDim == 2) {  ForceProj_Vector[0] += Weight_ObjFunc*(y - y_origin)/RefLength; ForceProj_Vector[1] += -Weight_ObjFunc*(x - x_origin)/RefLength; }
            if (nDim == 3) { ForceProj_Vector[0] += Weight_ObjFunc*(y - y_origin)/RefLength; ForceProj_Vector[1] += -Weight_ObjFunc*(x - x_origin)/RefLength; ForceProj_Vector[2] += 0; }
            break;
          case FIGURE_OF_MERIT :
            if (nDim == 2) { SU2_MPI::Error("This functional is not possible in 2D!!", CURRENT_FUNCTION);}
            if (nDim == 3) {
              ForceProj_Vector[0] += -Weight_ObjFunc*invCQ;
              ForceProj_Vector[1] += -Weight_ObjFunc*CTRCQ2*(z - z_origin);
              ForceProj_Vector[2] +=  Weight_ObjFunc*CTRCQ2*(y - y_origin);
            }
            break;

          default :
            break;
        }
        
        /*--- Store the force projection vector at this node ---*/
        
        nodes->SetForceProj_Vector(iPoint,ForceProj_Vector);
        
      }
    }
  }
  delete [] ForceProj_Vector;
  
}

void CAdjEulerSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  unsigned short iMarker, iVar, jVar, kVar, iDim, jDim, iIndex;
  unsigned long iVertex, iPoint, iPointNearField, nPointNearField = 0;
  su2double factor = 1.0, AngleDouble, data, aux, *IntBound_Vector, *coord, *FlowSolution, WeightSB, MinDist = 1E6, Dist, DerivativeOF = 0.0, *Normal, Area, UnitNormal[3], velocity[3], Energy, Rho, sqvel, proj_vel, phi, a1, a2;
  su2double **A, **M, **AM, *b;
  short AngleInt = 0, IndexNF_inv[180], iColumn;
  ifstream index_file;
  string text_line;
  vector<vector<su2double> > NearFieldWeight;
  vector<su2double> CoordNF;
  vector<short> IndexNF;
  
  IntBound_Vector = new su2double [nVar];
  
  /*--- Allocate vectors and matrices ---*/
  
  b = new su2double [nVar];
  A = new su2double* [nVar];
  M = new su2double* [nVar];
  AM = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    A[iVar] = new su2double [nVar];
    M[iVar] = new su2double [nVar];
    AM[iVar] = new su2double [nVar];
  }
  
  /*--- If equivalent area objective function, read the value of
   the derivative from a file, this is a preprocess of the direct solution ---*/
  
  if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) {
    
    /*--- Read derivative of the objective function at the NearField from file ---*/
    index_file.open("WeightNF.dat", ios::in);
    if (index_file.fail()) {
      SU2_MPI::Error("There is no Weight Nearfield Pressure file (WeightNF.dat).", CURRENT_FUNCTION);
    }
    
    nPointNearField = 0;
    
    while (index_file) {
      string line;
      getline(index_file, line);
      istringstream is(line);
      
      /*--- The first row provides the azimuthal angle ---*/
      
      if (nPointNearField == 0) {
        is >> data; // The first column is related with the coordinate
        while (is.good()) { is >> data; IndexNF.push_back(SU2_TYPE::Int(data)); }
      }
      else {
        is >> data; CoordNF.push_back(data); // The first column is the point coordinate
        vector<su2double> row;
        while (is.good()) { is >> data; row.push_back(data); }
        NearFieldWeight.push_back(row);
      }
      nPointNearField++;
    }
    
    /*--- Note tha the first row is the azimuthal angle ---*/
    
    nPointNearField = nPointNearField - 1;
    
    for (AngleInt = 0; AngleInt < 180; AngleInt++)
      IndexNF_inv[AngleInt] = -1;
    
  if (IndexNF.size() <= 180) {
    for (iIndex = 0; iIndex < IndexNF.size(); iIndex++)
      IndexNF_inv[IndexNF[iIndex]] = iIndex;
  }
  else {
    SU2_MPI::Error("", CURRENT_FUNCTION);
  }
    
  }
  
  /*--- Compute the jump on the adjoint variables for the upper and the lower side ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);
        
        for (iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = Normal[iDim]/Area;
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          coord = geometry->node[iPoint]->GetCoord();
          DerivativeOF = 0.0;
          
          /*--- Just in case the functional depend also on the surface pressure ---*/
          
          WeightSB = 1.0-config->GetWeightCd();
          
          su2double AoA, XcoordRot = 0.0, YcoordRot = 0.0, ZcoordRot = 0.0;
          
          if (nDim == 2) XcoordRot = coord[0];
          if (nDim == 3) {
            
            /*--- Rotate the nearfield cylinder  ---*/
            
            AoA = -(config->GetAoA()*PI_NUMBER/180.0);
            XcoordRot = coord[0]*cos(AoA) - coord[2]*sin(AoA);
            YcoordRot = coord[1];
            ZcoordRot = coord[0]*sin(AoA) + coord[2]*cos(AoA);
          }
          
          switch (config->GetKind_ObjFunc()) {
            case EQUIVALENT_AREA :
              
              if (nDim == 2) AngleInt = 0;
              
              if (nDim == 3) {
                
                /*--- Compute the azimuthal angle of the iPoint ---*/
                
                AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);
                
                /*--- Fix an azimuthal line due to misalignments of the near-field ---*/
                
                su2double FixAzimuthalLine = config->GetFixAzimuthalLine();
                
                if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1)) AngleDouble = FixAzimuthalLine - 0.1;
                
                AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));
                if (AngleInt < 0) AngleInt = 180 + AngleInt;
                
              }
              
              if (AngleInt <= 60) {
                iColumn = IndexNF_inv[AngleInt];
                
                /*--- An azimuthal angle is not defined... this happens with MG levels ---*/
                
                if (iColumn < 0.0) {
                  if (IndexNF_inv[AngleInt+1] > 0) { iColumn = IndexNF_inv[AngleInt+1]; goto end; }
                  if (IndexNF_inv[AngleInt-1] > 0) { iColumn = IndexNF_inv[AngleInt-1]; goto end; }
                  if (IndexNF_inv[AngleInt+2] > 0) { iColumn = IndexNF_inv[AngleInt+2]; goto end; }
                  if (IndexNF_inv[AngleInt-2] > 0) { iColumn = IndexNF_inv[AngleInt-2]; goto end; }
                  if (IndexNF_inv[AngleInt+3] > 0) { iColumn = IndexNF_inv[AngleInt+3]; goto end; }
                  if (IndexNF_inv[AngleInt-3] > 0) { iColumn = IndexNF_inv[AngleInt-3]; goto end; }
                  if (IndexNF_inv[AngleInt+4] > 0) { iColumn = IndexNF_inv[AngleInt+4]; goto end; }
                  if (IndexNF_inv[AngleInt-4] > 0) { iColumn = IndexNF_inv[AngleInt-4]; goto end; }
                }
                
              end:
                
                if (iColumn < 0.0) { cout <<" An azimuthal angle is not defined..." << endl; }
                
                /*--- Find the value of the weight in the table, using the azimuthal angle  ---*/
                
                MinDist = 1E6;
                for (iPointNearField = 0; iPointNearField < nPointNearField; iPointNearField++) {
                  Dist = fabs(CoordNF[iPointNearField] - XcoordRot);
                  if (Dist <= MinDist) {
                    MinDist = Dist;
                    DerivativeOF = factor*WeightSB*NearFieldWeight[iPointNearField][iColumn];
                  }
                }
              }
              else DerivativeOF = 0.0;
              
              if ((MinDist > 1E-6) || (coord[nDim-1] > 0.0)) DerivativeOF = 0.0;
              
              break;
              
            case NEARFIELD_PRESSURE :
              
              DerivativeOF = factor*WeightSB*(solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint)
                                              - solver_container[FLOW_SOL]->GetPressure_Inf());
              
              break;
              
          }
          
          /*--- Compute the jump of the adjoint variables (2D, and 3D problems) --*/
          
          FlowSolution = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          
          Rho = FlowSolution[0];
          Energy = FlowSolution[nVar-1]/FlowSolution[0];
          
          sqvel = 0.0; proj_vel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            velocity[iDim] = FlowSolution[iDim+1]/FlowSolution[0];
            sqvel    += velocity[iDim]*velocity[iDim];
            proj_vel += velocity[iDim]*UnitNormal[iDim];
          }
          
          if (nDim == 2) {
            
            /*--- Compute the projected Jacobian ---*/
            
            A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 1.0; A[0][3] = 0.0;
            A[1][0] = -velocity[0]*velocity[1]; A[1][1] = velocity[1]; A[1][2] = velocity[0]; A[1][3] = 0.0;
            A[2][0] = 0.5*(Gamma-3.0)*velocity[1]*velocity[1]+0.5*Gamma_Minus_One*velocity[0]*velocity[0]; A[2][1] = -Gamma_Minus_One*velocity[0];
            A[2][2] = (3.0-Gamma)*velocity[1]; A[2][3] = Gamma_Minus_One; A[3][0] = -Gamma*velocity[1]*Energy+Gamma_Minus_One*velocity[1]*sqvel;
            A[3][1] = -Gamma_Minus_One*velocity[0]*velocity[1]; A[3][2] = Gamma*Energy-0.5*Gamma_Minus_One*(velocity[0]*velocity[0]+3.0*velocity[1]*velocity[1]);  A[3][3] = Gamma*velocity[1];
            
            /*--- Compute the transformation matrix ---*/
            
            M[0][0] = 1.0; M[0][1] = 0.0; M[0][2] = 0.0; M[0][3] = 0.0;
            M[1][0] = velocity[0]; M[1][1] = Rho; M[1][2] = 0.0; M[1][3] = 0.0;
            M[2][0] = velocity[1]; M[2][1] = 0.0; M[2][2] = Rho; M[2][3] = 0.0;
            M[3][0] = 0.5*sqvel;  M[3][1] = Rho*velocity[0]; M[3][2] = Rho*velocity[1]; M[3][3] = 1.0/Gamma_Minus_One;
            
            /*--- Create the soruce term (AM)^T X = b ---*/
            
            b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = DerivativeOF;
            
          }
          
          if (nDim == 3) {
            
            
            /*--- Compute the projected Jacobian ---*/
            
            phi = 0.5*Gamma_Minus_One*sqvel;
            a1 = Gamma*Energy-phi; a2 = Gamma-1.0;
            
            A[0][0] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) A[0][iDim+1] = UnitNormal[iDim];
            A[0][nDim+1] = 0.0;
            
            for (iDim = 0; iDim < nDim; iDim++) {
              A[iDim+1][0] = (UnitNormal[iDim]*phi - velocity[iDim]*proj_vel);
              for (jDim = 0; jDim < nDim; jDim++)
                A[iDim+1][jDim+1] = (UnitNormal[jDim]*velocity[iDim]-a2*UnitNormal[iDim]*velocity[jDim]);
              A[iDim+1][iDim+1] += proj_vel;
              A[iDim+1][nDim+1] = a2*UnitNormal[iDim];
            }
            
            A[nDim+1][0] = proj_vel*(phi-a1);
            for (iDim = 0; iDim < nDim; iDim++)
              A[nDim+1][iDim+1] = (UnitNormal[iDim]*a1-a2*velocity[iDim]*proj_vel);
            A[nDim+1][nDim+1] = Gamma*proj_vel;
            
            /*--- Compute the transformation matrix ---*/
            
            M[0][0] = 1.0; M[0][1] = 0.0; M[0][2] = 0.0; M[0][3] = 0.0; M[0][4] = 0.0;
            M[1][0] = velocity[0]; M[1][1] = Rho; M[1][2] = 0.0; M[1][3] = 0.0; M[1][4] = 0.0;
            M[2][0] = velocity[1]; M[2][1] = 0.0; M[2][2] = Rho; M[2][3] = 0.0; M[2][4] = 0.0;
            M[3][0] = velocity[2]; M[3][1] = 0.0; M[3][2] = 0.0; M[3][3] = Rho; M[3][4] = 0.0;
            M[4][0] = 0.5*sqvel; M[4][1] = Rho*velocity[0]; M[4][2] = Rho*velocity[1];
            M[4][3] = Rho*velocity[2]; M[4][4] = 1.0/Gamma_Minus_One;
            
            /*--- Create the soruce term (AM)^T X = b ---*/
            
            b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; b[4] = DerivativeOF;
            
          }
          
          /*--- Compute A times M ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++) {
              aux = 0.0;
              for (kVar = 0; kVar < nVar; kVar++)
                aux += A[iVar][kVar]*M[kVar][jVar];
              AM[iVar][jVar] = aux;
            }
          
          /*--- Compute the transpose matrix ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              A[iVar][jVar] = AM[jVar][iVar];
          
          /*--- Solve the linear system using a LU descomposition --*/
          
          Gauss_Elimination(A, b, nVar);
          
          /*--- Update the internal boundary jump --*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            IntBound_Vector[iVar] = b[iVar];
          
          nodes->SetIntBoundary_Jump(iPoint,IntBound_Vector);
          
        }
      }
  
  delete [] IntBound_Vector;
  
  /*--- Deallocate the linear system ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] A[iVar];
    delete [] M[iVar];
    delete [] AM[iVar];
  }
  delete [] A;
  delete [] M;
  delete [] AM;
  delete [] b;
  
}

void CAdjEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution, *Solution_Fine;
  
  bool restart = config->GetRestart();
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  if (restart) {
    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][ADJFLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][ADJFLOW_SOL]->GetNodes()->SetSolution(iPoint, Solution);
        
      }
      solver_container[iMesh][ADJFLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][ADJFLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    }
    delete [] Solution;
  }
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if ((TimeIter == 0) && (dual_time)) {
      solver_container[iMesh][ADJFLOW_SOL]->GetNodes()->Set_Solution_time_n();
      solver_container[iMesh][ADJFLOW_SOL]->GetNodes()->Set_Solution_time_n1();
    }
  }
  
}

void CAdjEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double SharpEdge_Distance;
  bool RightSol = true;
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool muscl          = config->GetMUSCL_AdjFlow();
  bool limiter        = (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);
  bool fixed_cl       = config->GetFixed_CL_Mode();
  bool eval_dof_dcx   = config->GetEval_dOF_dCX();

  /*--- Update the objective function coefficient to guarantee zero gradient. ---*/
  
  if (fixed_cl && eval_dof_dcx) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Residual initialization ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Get the distance form a sharp edge ---*/
    
    SharpEdge_Distance = geometry->node[iPoint]->GetSharpEdge_Distance();
    
    /*--- Initialize the non-physical points vector ---*/

    nodes->SetNon_Physical(iPoint,false);
    
    /*--- Set the primitive variables compressible
     adjoint variables ---*/
    
    RightSol = nodes->SetPrimVar(iPoint,SharpEdge_Distance, false, config);
    if (!RightSol) { nodes->SetNon_Physical(iPoint,true); ErrorCounter++; }
    
    /*--- Initialize the convective residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  
  if ((muscl) && (iMesh == MESH_0)) {
    
    /*--- Compute gradients for upwind second-order reconstruction ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
    
    if (limiter && !Output) SetSolution_Limiter(geometry, config);
    
  }
  
  /*--- Artificial dissipation for centered schemes ---*/
  
  if (center) {
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    }
  }
  
  /*--- Initialize the Jacobian for implicit integration ---*/
  
  if (implicit) Jacobian.SetValZero();
  
  /*--- Error message ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CAdjEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool jst_scheme = ((config->GetKind_Centered_AdjFlow() == JST) && (iMesh == MESH_0));
  bool grid_movement  = config->GetGrid_Movement();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, normal, and neighbors---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Adjoint variables w/o reconstruction ---*/

    numerics->SetAdjointVar(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    
    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetConservative(solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint),
                              solver_container[FLOW_SOL]->GetNodes()->GetSolution(jPoint));
    
    numerics->SetSoundSpeed(solver_container[FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint),
                            solver_container[FLOW_SOL]->GetNodes()->GetSoundSpeed(jPoint));
    numerics->SetEnthalpy(solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint),
                          solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(jPoint));

    
    numerics->SetLambda(solver_container[FLOW_SOL]->GetNodes()->GetLambda(iPoint),
                        solver_container[FLOW_SOL]->GetNodes()->GetLambda(jPoint));
    
    if (jst_scheme) {
      numerics->SetUndivided_Laplacian(nodes->GetUndivided_Laplacian(iPoint), nodes->GetUndivided_Laplacian(jPoint));
      numerics->SetSensor(nodes->GetSensor(iPoint), nodes->GetSensor(jPoint));
    }
    
    /*--- Mesh motion ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals ---*/

    numerics->ComputeResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j,
                              Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/

    LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
    LinSysRes.SubtractBlock(jPoint, Res_Conv_j);
    LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
    LinSysRes.SubtractBlock(jPoint, Res_Visc_j);
    
    /*--- Implicit contribution to the residual ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
}


void CAdjEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL,
  *Limiter_j = NULL, *Psi_i = NULL, *Psi_j = NULL, *V_i, *V_j, Non_Physical = 1.0;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool implicit         = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_AdjFlow() && (iMesh == MESH_0));
  bool limiter          = (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER);
  bool grid_movement    = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Adjoint variables w/o reconstruction ---*/
    
    Psi_i = nodes->GetSolution(iPoint);
    Psi_j = nodes->GetSolution(jPoint);
    numerics->SetAdjointVar(Psi_i, Psi_j);
    
    /*--- Primitive variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
    V_j = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(jPoint);
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Grid velocities for dynamic meshes ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (muscl) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Adjoint variables using gradient reconstruction and limiters ---*/

      Gradient_i = nodes->GetGradient(iPoint); Gradient_j = nodes->GetGradient(jPoint);
      if (limiter) { Limiter_i = nodes->GetLimiter(iPoint); Limiter_j = nodes->GetLimiter(jPoint); }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0; Project_Grad_j = 0;
        Non_Physical = nodes->GetNon_Physical(iPoint)*nodes->GetNon_Physical(jPoint);
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
        }
        if (limiter) {
          Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i*Limiter_i[iDim];
          Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j*Limiter_j[iDim];
        }
        else {
          Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetAdjointVar(Solution_i, Solution_j);
      
    }
    
    /*--- Compute the residual---*/
    
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Add and Subtract Residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual_i);
    LinSysRes.SubtractBlock(jPoint, Residual_j);
    
    /*--- Implicit contribution to the residual ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
  
}

void CAdjEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                      CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  bool axisymmetric   = config->GetAxisymmetric();
  //  bool gravity        = (config->GetGravityForce() == YES);
  bool harmonic_balance  = (config->GetTime_Marching() == HARMONIC_BALANCE);
  
  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
  if (rotating_frame) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the adjoint variables ---*/
      numerics->SetAdjointVar(nodes->GetSolution(iPoint),
                              nodes->GetSolution(iPoint));
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the adjoint rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  if (harmonic_balance) {
    
    su2double Volume, Source;
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Get stored harmonic balance source term ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Source = nodes->GetHarmonicBalance_Source(iPoint,iVar);
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution ---*/
      numerics->SetConservative(solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint), solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint));
      
      /*--- Set adjoint variables ---*/
      numerics->SetAdjointVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));
      
      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set coordinate ---*/
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  //  if (gravity) {
  //
  //  }
  
}

void CAdjEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                      CConfig *config, unsigned short iMesh) {
}

void CAdjEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge;
  unsigned short iVar;
  su2double *Diff;
  
  Diff = new su2double[nVar];
  
  nodes->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = nodes->GetSolution(iPoint,iVar) - nodes->GetSolution(jPoint,iVar);
    
#ifdef STRUCTURED_GRID
    
    if (geometry->node[iPoint]->GetDomain()) nodes->SubtractUnd_Lapl(iPoint, Diff);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddUnd_Lapl(jPoint, Diff);
    
#else
    
    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) nodes->SubtractUnd_Lapl(iPoint, Diff);
      if (geometry->node[jPoint]->GetDomain()) nodes->AddUnd_Lapl(jPoint, Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) nodes->SubtractUnd_Lapl(iPoint, Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) nodes->AddUnd_Lapl(jPoint, Diff);
    
#endif
    
  }
  
#ifdef STRUCTURED_GRID
  
  unsigned long Point_Normal = 0, iVertex;
  unsigned short iMarker;
  su2double *Psi_mirror;
  
  Psi_mirror = new su2double[nVar];
  
  /*--- Loop over all boundaries and include an extra contribution
   from a halo node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
          
          /*--- Interpolate & compute difference in the conserved variables ---*/
          
          for (iVar = 0; iVar < nVar; iVar++) {
            Psi_mirror[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - nodes->GetSolution(Point_Normal, iVar);
            Diff[iVar]   = nodes->GetSolution(iPoint,iVar) - Psi_mirror[iVar];
          }
          
          /*--- Subtract contribution at the boundary node only ---*/
          
          nodes->SubtractUnd_Lapl(iPoint,Diff);
        }
      }
    }
  }
  
  delete [] Psi_mirror;
  
#endif
  
  delete [] Diff;
  
  /*--- MPI parallelization ---*/
  
  InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
  CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);
  
}

void CAdjEulerSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint;
  su2double SharpEdge_Distance, eps, ds, scale, Sensor, Param_Kappa_2, Param_Kappa_4;
  
  eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
  Param_Kappa_2 = config->GetKappa_2nd_AdjFlow();
  Param_Kappa_4 = config->GetKappa_4th_AdjFlow();
  
  if (Param_Kappa_2 != 0.0) scale = 2.0 * Param_Kappa_4 / Param_Kappa_2;
  else scale = 0.0;
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    SharpEdge_Distance = (geometry->node[iPoint]->GetSharpEdge_Distance() - config->GetAdjSharp_LimiterCoeff()*eps);
    
    ds = 0.0;
    if (SharpEdge_Distance < -eps) ds = 1.0;
    if (fabs(SharpEdge_Distance) <= eps) ds = 1.0 - (0.5*(1.0+(SharpEdge_Distance/eps)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps)));
    if (SharpEdge_Distance > eps) ds = 0.0;
    
    Sensor = scale * ds;
    
    nodes->SetSensor(iPoint,Sensor);
    
  }
  
  /*--- MPI parallelization ---*/
  
  InitiateComms(geometry, config, SENSOR);
  CompleteComms(geometry, config, SENSOR);
  
}

void CAdjEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                           CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) / Vol;
    
    Res_TruncError = nodes->GetResTruncError(iPoint);
    Residual = LinSysRes.GetBlock(iPoint);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      nodes->AddSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
    
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) / Vol;
    
    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    local_Residual = LinSysRes.GetBlock(iPoint);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Res = local_Residual[iVar] + local_Res_TruncError[iVar];
      nodes->AddSolution(iPoint,iVar, -Res*Delta);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
    
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, *local_Res_TruncError, Vol;
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    if (solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) != 0.0) {
      Delta = Vol / solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint);
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = -(LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
    
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    for (iVar = 0; iVar < nVar; iVar++) {
      nodes->AddSolution(iPoint,iVar, config->GetRelaxation_Factor_AdjFlow()*LinSysSol[iPoint*nVar+iVar]);
    }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CAdjEulerSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  unsigned long iVertex, iPoint, Neigh;
  unsigned short iPos, jPos;
  unsigned short iDim, iMarker, iNeigh;
  su2double *d = NULL, *Normal = NULL, *Psi = NULL, *U = NULL, Enthalpy, conspsi = 0.0, Mach_Inf,
  Area, **PrimVar_Grad = NULL, *ConsPsi_Grad = NULL,
  ConsPsi, d_press, grad_v, v_gradconspsi, UnitNormal[3], *GridVel = NULL,
  eps, r, ru, rv, rw, rE, p, T, dp_dr, dp_dru, dp_drv,
  dp_drw, dp_drE, dH_dr, dH_dru, dH_drv, dH_drw, dH_drE, H, *USens, D[3][3], Dd[3], scale = 1.0;
  su2double RefVel2, RefDensity, Mach2Vel, *Velocity_Inf, factor;
  su2double Vn, SoundSpeed, *Velocity;
  
  USens = new su2double[nVar];
  Velocity = new su2double[nDim];

  su2double Gas_Constant    = config->GetGas_ConstantND();
  bool grid_movement     = config->GetGrid_Movement();
  su2double RefArea    = config->GetRefArea();
  su2double Mach_Motion     = config->GetMach_Motion();
  unsigned short ObjFunc = config->GetKind_ObjFunc();

  if (config->GetSystemMeasurements() == US) scale = 1.0/12.0;
  else scale = 1.0;
  
  /*--- Compute non-dimensional factor. For dynamic meshes, use the motion Mach 
   number as a reference value for computing the force coefficients. 
   Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  RefDensity  = config->GetDensity_FreeStreamND();

  factor = 1.0;
  /*-- For multi-objective problems these scaling factors are applied before solution ---*/
  if (config->GetnObj()==1) {
    factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);
    if ((ObjFunc == INVERSE_DESIGN_HEATFLUX)  ||
        (ObjFunc == TOTAL_HEATFLUX) || (ObjFunc == MAXIMUM_HEATFLUX) ||
  (ObjFunc == SURFACE_MASSFLOW) )
      factor = 1.0;

    if ((ObjFunc == SURFACE_TOTAL_PRESSURE) || (ObjFunc == SURFACE_STATIC_PRESSURE)) factor = 1.0/Area_Monitored;

  }

  /*--- Initialize sensitivities to zero ---*/
  
  Total_Sens_Geo = 0.0;
  Total_Sens_Mach = 0.0;
  Total_Sens_AoA = 0.0;
  Total_Sens_Press = 0.0;
  Total_Sens_Temp = 0.0;
  Total_Sens_BPress = 0.0;
  
  /*--- Loop over boundary markers to select those for Euler walls ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)

    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)
      
    /*--- Loop over points on the surface to store the auxiliary variable ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);
          conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
          for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
          
          nodes->SetAuxVar(iPoint,conspsi);
          
          /*--- Also load the auxiliary variable for first neighbors ---*/
          
          for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
            Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
            Psi = nodes->GetSolution(Neigh);
            U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(Neigh);
            Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(Neigh);
            conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
            for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
            nodes->SetAuxVar(Neigh,conspsi);
          }
        }
      }
  
  /*--- Compute surface gradients of the auxiliary variable ---*/
  
  SetAuxVar_Surface_Gradient(geometry, config);
  
  /*--- Evaluate the shape sensitivity ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker] = 0.0;
    
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          
          d = nodes->GetForceProj_Vector(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          
          PrimVar_Grad = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
          ConsPsi_Grad = nodes->GetAuxVarGradient(iPoint);
          ConsPsi = nodes->GetAuxVar(iPoint);
          
          d_press = 0.0; grad_v = 0.0; v_gradconspsi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*-- Retrieve the value of the pressure gradient ---*/
            
            d_press += d[iDim]*PrimVar_Grad[nDim+1][iDim];
            
            /*-- Retrieve the value of the velocity gradient ---*/
            
            grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
            
            /*-- Retrieve the value of the theta gradient ---*/
            
            v_gradconspsi += solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim) * ConsPsi_Grad[iDim];
            
            /*--- Additional sensitivity term for grid movement ---*/
            
            if (grid_movement) {
              GridVel = geometry->node[iPoint]->GetGridVel();
              v_gradconspsi -= GridVel[iDim] * ConsPsi_Grad[iDim];
            }
            
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          
          CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area * scale * factor;
          
          /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
          
          if (config->GetSens_Remove_Sharp()) {
            eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
            if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
              CSensitivity[iMarker][iVertex] = 0.0;
          }
          
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex];
          
        }
      }
      
      Total_Sens_Geo += Sens_Geo[iMarker];
      
    }
  }
  
  
  /*--- Farfield Sensitivity (Mach, AoA, Press, Temp), only for compressible flows ---*/
  

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_BPress[iMarker] = 0.0;
      if (config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);

          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

          Vn = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = U[iDim+1]/U[0];
            Vn += UnitNormal[iDim]*Velocity[iDim];
          }

          SoundSpeed = solver_container[FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint);
            if (Vn<SoundSpeed && Vn>0) {
              /*TODO: MDO compatible*/
            Sens_BPress[iMarker]+=Psi[nDim+1]*(SoundSpeed*SoundSpeed-Vn*Vn)/(Vn*Gamma_Minus_One);
            if (config->GetKind_ObjFunc()==SURFACE_STATIC_PRESSURE)
              Sens_BPress[iMarker]+=1;
            if (config->GetKind_ObjFunc()==SURFACE_TOTAL_PRESSURE) {
              for (iDim=0; iDim<nDim; iDim++)
                Sens_BPress[iMarker]+=0.5*Velocity[iDim]*Velocity[iDim]/(Vn*Vn);
            }
          }
        }
        Total_Sens_BPress+= Sens_BPress[iMarker] * scale * factor;
      }
    }

    if (config->GetMarker_All_KindBC(iMarker) == FAR_FIELD || config->GetMarker_All_KindBC(iMarker) == INLET_FLOW
        || config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_INLET || config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_OUTLET
        || config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW  ) {

      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          r = U[0]; ru = U[1]; rv = U[2];
          if (nDim == 2) { rw = 0.0; rE = U[3]; }
          else { rw = U[3]; rE = U[4]; }
          p = Gamma_Minus_One*(rE-(ru*ru + rv*rv + rw*rw)/(2*r));

          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

          H = (rE + p)/r;

          dp_dr = Gamma_Minus_One*(ru*ru + rv*rv + rw*rw)/(2*r*r);
          dp_dru = -Gamma_Minus_One*ru/r;
          dp_drv = -Gamma_Minus_One*rv/r;
          if (nDim == 2) { dp_drw = 0.0; dp_drE = Gamma_Minus_One; }
          else { dp_drw = -Gamma_Minus_One*rw/r; dp_drE = Gamma_Minus_One; }

          dH_dr = (-H + dp_dr)/r; dH_dru = dp_dru/r; dH_drv = dp_drv/r;
          if (nDim == 2) { dH_drw = 0.0; dH_drE = (1 + dp_drE)/r; }
          else { dH_drw = dp_drw/r; dH_drE = (1 + dp_drE)/r; }

          if (nDim == 2) {
            Jacobian_j[0][0] = 0.0;
            Jacobian_j[1][0] = Area*UnitNormal[0];
            Jacobian_j[2][0] = Area*UnitNormal[1];
            Jacobian_j[3][0] = 0.0;

            Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitNormal[0] + (-(ru*rv)/(r*r))*Area*UnitNormal[1];
            Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1];
            Jacobian_j[2][1] = (dp_drv)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[1];
            Jacobian_j[3][1] = (dp_drE)*Area*UnitNormal[0];

            Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitNormal[0] + (-(rv*rv)/(r*r) + dp_dr)*Area*UnitNormal[1];
            Jacobian_j[1][2] = (rv/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[1];
            Jacobian_j[2][2] = (ru/r)*Area*UnitNormal[0] + (2*rv/r + dp_drv)*Area*UnitNormal[1];
            Jacobian_j[3][2] = (dp_drE)*Area*UnitNormal[1];

            Jacobian_j[0][3] = (ru*dH_dr)*Area*UnitNormal[0] + (rv*dH_dr)*Area*UnitNormal[1];
            Jacobian_j[1][3] = (H + ru*dH_dru)*Area*UnitNormal[0] + (rv*dH_dru)*Area*UnitNormal[1];
            Jacobian_j[2][3] = (ru*dH_drv)*Area*UnitNormal[0] + (H + rv*dH_drv)*Area*UnitNormal[1];
            Jacobian_j[3][3] = (ru*dH_drE)*Area*UnitNormal[0] + (rv*dH_drE)*Area*UnitNormal[1];
          }
          else {
            Jacobian_j[0][0] = 0.0;
            Jacobian_j[1][0] = Area*UnitNormal[0];
            Jacobian_j[2][0] = Area*UnitNormal[1];
            Jacobian_j[3][0] = Area*UnitNormal[2];
            Jacobian_j[4][0] = 0.0;

            Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitNormal[0] + (-(ru*rv)/(r*r))*Area*UnitNormal[1] + (-(ru*rw)/(r*r))*Area*UnitNormal[2];
            Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1] + (rw/r)*Area*UnitNormal[2];
            Jacobian_j[2][1] = (dp_drv)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[1];
            Jacobian_j[3][1] = (dp_drw)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[2];
            Jacobian_j[4][1] = (dp_drE)*Area*UnitNormal[0];

            Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitNormal[0] + (-(rv*rv)/(r*r) + dp_dr)*Area*UnitNormal[1] + (-(rv*rw)/(r*r))*Area*UnitNormal[2];
            Jacobian_j[1][2] = (rv/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[1];
            Jacobian_j[2][2] = (ru/r)*Area*UnitNormal[0] + (2*rv/r + dp_drv)*Area*UnitNormal[1] + (rw/r)*Area*UnitNormal[2];
            Jacobian_j[3][2] = (dp_drw)*Area*UnitNormal[1] + (rv/r)*Area*UnitNormal[2];
            Jacobian_j[4][2] = (dp_drE)*Area*UnitNormal[1];

            Jacobian_j[0][3] = (-(ru*rw)/(r*r))*Area*UnitNormal[0] + (-(rv*rw)/(r*r))*Area*UnitNormal[1] + (-(rw*rw)/(r*r) + dp_dr)*Area*UnitNormal[2];
            Jacobian_j[1][3] = (rw/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[2];
            Jacobian_j[2][3] = (rw/r)*Area*UnitNormal[1] + (dp_drv)*Area*UnitNormal[2];
            Jacobian_j[3][3] = (ru/r)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1] + (2*rw/r + dp_drw)*Area*UnitNormal[2];
            Jacobian_j[4][3] = (dp_drE)*Area*UnitNormal[2];

            Jacobian_j[0][4] = (ru*dH_dr)*Area*UnitNormal[0] + (rv*dH_dr)*Area*UnitNormal[1] + (rw*dH_dr)*Area*UnitNormal[2];
            Jacobian_j[1][4] = (H + ru*dH_dru)*Area*UnitNormal[0] + (rv*dH_dru)*Area*UnitNormal[1] + (rw*dH_dru)*Area*UnitNormal[2];
            Jacobian_j[2][4] = (ru*dH_drv)*Area*UnitNormal[0] + (H + rv*dH_drv)*Area*UnitNormal[1] + (rw*dH_drv)*Area*UnitNormal[2];
            Jacobian_j[3][4] = (ru*dH_drw)*Area*UnitNormal[0] + (rv*dH_drw)*Area*UnitNormal[1] + (H + rw*dH_drw)*Area*UnitNormal[2];
            Jacobian_j[4][4] = (ru*dH_drE)*Area*UnitNormal[0] + (rv*dH_drE)*Area*UnitNormal[1] + (rw*dH_drE)*Area*UnitNormal[2];
          }

          /*--- Mach number sensitivity ---*/

          USens[0] = 0.0; USens[1] = ru/Mach_Inf; USens[2] = rv/Mach_Inf;
          if (nDim == 2) { USens[3] = Gamma*Mach_Inf*p; }
          else { USens[3] = rw/Mach_Inf; USens[4] = Gamma*Mach_Inf*p; }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_Mach[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }

          /*--- AoA sensitivity ---*/

          USens[0] = 0.0;
          if (nDim == 2) { USens[1] = -rv; USens[2] = ru; USens[3] = 0.0; }
          else { USens[1] = -rw; USens[2] = 0.0; USens[3] = ru; USens[4] = 0.0; }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_AoA[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }

          /*--- Pressure sensitivity ---*/

          USens[0] = r/p; USens[1] = ru/p; USens[2] = rv/p;
          if (nDim == 2) { USens[3] = rE/p; }
          else { USens[3] = rw/p; USens[4] = rE/p; }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_Press[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }

          /*--- Temperature sensitivity ---*/

          T = p/(r*Gas_Constant);
          USens[0] = -r/T; USens[1] = 0.5*ru/T; USens[2] = 0.5*rv/T;
          if (nDim == 2) { USens[3] = (ru*ru + rv*rv + rw*rw)/(r*T); }
          else { USens[3] = 0.5*rw/T; USens[4] = (ru*ru + rv*rv + rw*rw)/(r*T); }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_Temp[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }
        }
      }

      Total_Sens_Mach  -= Sens_Mach[iMarker] * scale * factor;
      Total_Sens_AoA   -= Sens_AoA[iMarker] * scale * factor;
      Total_Sens_Press -= Sens_Press[iMarker] * scale * factor;
      Total_Sens_Temp  -= Sens_Temp[iMarker] * scale * factor;

    }
  }

  /*--- Explicit contribution from objective function quantity ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {

      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          p = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          d = nodes->GetForceProj_Vector(iPoint);
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

          /*--- Mach number sensitivity ---*/

          for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = -(2.0/Mach_Inf)*d[iPos];
          for (iPos = 0; iPos < nDim; iPos++) Sens_Mach[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

          /*--- AoA sensitivity ---*/

          if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT ||
              config->GetKind_ObjFunc() == LIFT_COEFFICIENT ||
              config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT ||
              config->GetKind_ObjFunc() == EQUIVALENT_AREA ||
              config->GetKind_ObjFunc() == NEARFIELD_PRESSURE) {
            if (nDim == 2) {
              D[0][0] = 0.0; D[0][1] = -1.0;
              D[1][0] = 1.0; D[1][1] = 0.0;
            }
            else {
              D[0][0] = 0.0; D[0][1] = 0.0; D[0][2] = -1.0;
              D[1][0] = 0.0; D[1][1] = 0.0; D[1][2] = 0.0;
              D[2][0] = 1.0; D[2][1] = 0.0; D[2][2] = 0.0;
            }
            for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = 0.0;
            for (iPos = 0; iPos < nDim; iPos++) {
              for (jPos = 0; jPos < nDim; jPos++)
                Dd[iPos] += D[iPos][jPos]*d[jPos];
            }
          }

          /*--- Coefficients with no explicit AoA dependece ---*/

          else {
            for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
          }

          for (iPos = 0; iPos < nDim; iPos++)
            Sens_AoA[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

          /*--- Pressure sensitivity ---*/

          for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = -(1/p)*d[iPos];
          for (iPos = 0; iPos<nDim; iPos++)
            Sens_Press[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

          /*--- Temperature sensitivity ---*/

          for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
          for (iPos = 0; iPos<nDim; iPos++)
            Sens_Temp[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

        }
      }

      Total_Sens_Mach  += Sens_Mach[iMarker] * scale * factor;
      Total_Sens_AoA   += Sens_AoA[iMarker] * scale * factor;
      Total_Sens_Press += Sens_Press[iMarker] * scale * factor;
      Total_Sens_Temp  += Sens_Temp[iMarker] * scale * factor;

    }
  }

  
#ifdef HAVE_MPI
  
  su2double MyTotal_Sens_Geo   = Total_Sens_Geo;     Total_Sens_Geo = 0.0;
  su2double MyTotal_Sens_Mach  = Total_Sens_Mach;    Total_Sens_Mach = 0.0;
  su2double MyTotal_Sens_AoA   = Total_Sens_AoA;     Total_Sens_AoA = 0.0;
  su2double MyTotal_Sens_Press = Total_Sens_Press;   Total_Sens_Press = 0.0;
  su2double MyTotal_Sens_Temp  = Total_Sens_Temp;    Total_Sens_Temp = 0.0;
  su2double MyTotal_Sens_BPress  = Total_Sens_BPress;    Total_Sens_BPress = 0.0;
  
  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_BPress, &Total_Sens_BPress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif
  
  delete [] USens;
  delete [] Velocity;
  
}

void CAdjEulerSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  unsigned short iMarker;
  unsigned long iVertex, jVertex, nVertex, iPoint;
  su2double **A, *b, Sens, *ArchLength, *Coord_begin, *Coord_end, dist;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      nVertex = geometry->nVertex[iMarker];
      
      /*--- Allocate the linear system ---*/
      
      A = new su2double* [nVertex];
      b = new su2double [nVertex];
      ArchLength = new su2double [nVertex];
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        A[iVertex] = new su2double [nVertex];
      }
      
      /*--- Initialization ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        b[iVertex] = 0.0; ArchLength[iVertex] = 0.0;
        for (jVertex = 0; jVertex < nVertex; jVertex++)
          A[iVertex][jVertex] = 0.0;
      }
      
      /*--- Set the arch length ---*/
      
      ArchLength[0] = 0.0;
      for (iVertex = 1; iVertex < nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex-1]->GetNode();
        Coord_begin = geometry->node[iPoint]->GetCoord();
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Coord_end = geometry->node[iPoint]->GetCoord();
        dist = sqrt (pow( Coord_end[0]-Coord_begin[0], 2.0) + pow( Coord_end[1]-Coord_begin[1], 2.0));
        ArchLength[iVertex] = ArchLength[iVertex-1] + dist;
      }
      
      /*--- Remove the trailing edge effect ---*/
      
      su2double MinPosSens = 0.0; su2double MinNegSens = 0.0;
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        Sens = CSensitivity[iMarker][iVertex];
        if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.01) { MinNegSens = Sens; break; }
      }
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        Sens = CSensitivity[iMarker][iVertex];
        if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.99) { MinPosSens = Sens; break; }
      }
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        if (ArchLength[iVertex] < ArchLength[nVertex-1]*0.01)
          CSensitivity[iMarker][iVertex] = MinNegSens;
        if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.99)
          CSensitivity[iMarker][iVertex] = MinPosSens;
      }
      
      /*--- Set the right hand side of the system ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        b[iVertex] = CSensitivity[iMarker][iVertex];
      }
      
      /*--- Set the mass matrix ---*/
      
      su2double Coeff = 0.0, BackDiff = 0.0, ForwDiff = 0.0, CentDiff = 0.0;
      su2double epsilon = 5E-5;
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        
        if ((iVertex != nVertex-1) && (iVertex != 0)) {
          BackDiff = (ArchLength[iVertex]-ArchLength[iVertex-1]);
          ForwDiff = (ArchLength[iVertex+1]-ArchLength[iVertex]);
          CentDiff = (ArchLength[iVertex+1]-ArchLength[iVertex-1]);
        }
        if (iVertex == nVertex-1) {
          BackDiff = (ArchLength[nVertex-1]-ArchLength[nVertex-2]);
          ForwDiff = (ArchLength[0]-ArchLength[nVertex-1]);
          CentDiff = (ArchLength[0]-ArchLength[nVertex-2]);
        }
        if (iVertex == 0) {
          BackDiff = (ArchLength[0]-ArchLength[nVertex-1]);
          ForwDiff = (ArchLength[1]-ArchLength[0]);
          CentDiff = (ArchLength[1]-ArchLength[nVertex-1]);
        }
        
        Coeff = epsilon*2.0/(BackDiff*ForwDiff*CentDiff);
        
        A[iVertex][iVertex] = Coeff*CentDiff;
        
        if (iVertex != 0) A[iVertex][iVertex-1] = -Coeff*ForwDiff;
        else A[iVertex][nVertex-1] = -Coeff*ForwDiff;
        
        if (iVertex != nVertex-1) A[iVertex][iVertex+1] = -Coeff*BackDiff;
        else A[iVertex][0] = -Coeff*BackDiff;
        
      }
      
      /*--- Add the gradient value in the main diagonal ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++)
        A[iVertex][iVertex] += 1.0;
      
      /*--- Dirichlet boundary condition ---*/
      
      unsigned long iVertex = SU2_TYPE::Int(nVertex/2);
      A[iVertex][iVertex] = 1.0;
      A[iVertex][iVertex+1] = 0.0;
      A[iVertex][iVertex-1] = 0.0;
      
      Gauss_Elimination(A, b, (unsigned short)nVertex);
      
      /*--- Set the new value of the sensitiviy ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++)
        CSensitivity[iMarker][iVertex] = b[iVertex];
      
      /*--- Deallocate the linear system ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++)
        delete [] A[iVertex];
      delete [] A;
      delete [] b;
      delete [] ArchLength;
      
    }
  }
  
  
}

void CAdjEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                      CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned long Iter_Fixed_CL = config->GetUpdate_AoA_Iter_Limit();
  unsigned long InnerIter       = config->GetInnerIter();
  bool Update_AoA             = false;
  su2double dCL_dAlpha        = config->GetdCL_dAlpha()*180.0/PI_NUMBER;
  //unsigned long Update_Alpha  = config->GetUpdate_Alpha();
  
  //if (ExtIter == 0) AoA_Counter = 0;
  
  /*--- Only the fine mesh level should check the convergence criteria ---*/
  
  if ((iMesh == MESH_0) && Output) {
    
    /*--- Initialize the update flag to false ---*/
    
    Update_AoA = false;
    
    /*--- Reevaluate the lift derivative with respect to Angle of Attack
     at a fix number of iterations ---*/
    
    if ((InnerIter % Iter_Fixed_CL == 0) && (InnerIter != 0)) {
      //AoA_Counter++;
      //if ((AoA_Counter <= Update_Alpha)) Update_AoA = true;
      Update_AoA = true;
    }
    
    /*--- Store the update boolean for use on other mesh levels in the MG ---*/
    
    config->SetUpdate_AoA(Update_AoA);
    
  }
  
  else {
    Update_AoA = config->GetUpdate_AoA();
  }
  
  /*--- If we are within two digits of convergence in the CL coefficient,
   compute an updated value for the AoA at the farfield. We are iterating
   on the AoA in order to match the specified fixed lift coefficient. ---*/
  
  if (Update_AoA && Output) {
    
    /*--- Retrieve the old ACoeff ---*/
    
    if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) ACoeff_old = config->GetdCD_dCL();
    else if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) ACoeff_old = config->GetdCMx_dCL();
    else if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) ACoeff_old = config->GetdCMy_dCL();
    else if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) ACoeff_old = config->GetdCMz_dCL();
    else ACoeff_old = 0.0;

    /*--- Estimate the increment in the A coeff, (note that the slope is negative, a decrease in
     *   the CL derivative requires an increase in the A coeff ---*/
    
    /*--- A good estimation to d(dOF/dalpha)/dA_coeff is dCL_dAlpha ---*/
    
    ACoeff_inc =  (1.0/dCL_dAlpha)*Total_Sens_AoA;
    
    /*--- Compute a new value for the A coeff based on the fine mesh only (radians)---*/
    
    if (iMesh == MESH_0) { ACoeff = ACoeff_old + ACoeff_inc; }
    else {
      if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) ACoeff = config->GetdCD_dCL();
      else if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) ACoeff = config->GetdCMx_dCL();
      else if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) ACoeff = config->GetdCMy_dCL();
      else if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) ACoeff = config->GetdCMz_dCL();
      else ACoeff = 0.0;
    }
    
    /*--- Only the fine mesh stores the updated values for ACoeff in config ---*/
    
    if (iMesh == MESH_0) {
      if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT)  config->SetdCD_dCL(ACoeff);
      else if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) config->SetdCMx_dCL(ACoeff);
      else if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) config->SetdCMy_dCL(ACoeff);
      else if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) config->SetdCMz_dCL(ACoeff);
      else { config->SetdCD_dCL(0.0); config->SetdCMx_dCL(0.0); config->SetdCMy_dCL(0.0); config->SetdCMz_dCL(0.0); }
    }
    
    /*--- Compute the adjoint boundary condition ---*/
    
    SetForceProj_Vector(geometry, solver_container, config);
    
  }
  
  /*--- Output some information to the console with the headers ---*/
  
  bool write_heads = ((InnerIter % Iter_Fixed_CL == 0) && (InnerIter != 0));
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output) {
    cout.precision(7);
    cout.setf(ios::fixed, ios::floatfield);
    cout << endl << "-------------------------- Adjoint Fixed CL Mode -------------------------" << endl;
    cout << "Target dCL/dAlpha: 0.0 (1/deg), current dCL/dAlpha: " << Total_Sens_AoA*PI_NUMBER/180;
    if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT)  cout << " (1/deg), current dCD/dCL: " << config->GetdCD_dCL() <<" "<< endl;
    else if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) cout << " (1/deg), current dCMx/dCL: " << config->GetdCMx_dCL() <<" "<< endl;
    else if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) cout << " (1/deg), current dCMy/dCL: " << config->GetdCMy_dCL() <<" "<< endl;
    else if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) cout << " (1/deg), current dCMz/dCL: " << config->GetdCMz_dCL() <<" "<< endl;
    else cout << endl;
    cout << "-------------------------------------------------------------------------" << endl << endl;
  }
  
  
}


void CAdjEulerSolver::BC_Euler_Wall(CGeometry      *geometry,
                                    CSolver        **solver_container,
                                    CNumerics      *conv_numerics,
                                    CNumerics      *visc_numerics,
                                    CConfig        *config,
                                    unsigned short val_marker) {

  unsigned long iVertex, iPoint;
  su2double *d = NULL, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitNormal;
  su2double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2;
  unsigned short iDim, iVar, jDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  UnitNormal = new su2double[nDim];
  Velocity = new su2double[nDim];
  Psi      = new su2double[nVar];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      /*--- Create a copy of the adjoint solution ---*/
      Psi_Aux = nodes->GetSolution(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
      
      /*--- Flow solution ---*/
      U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
      
      /*--- Read the value of the objective function ---*/
      d = nodes->GetForceProj_Vector(iPoint);
      
      /*--- Normal vector computation ---*/
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

      for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = U[iDim+1] / U[0];

      Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);
      sq_vel   = 0.5*solver_container[FLOW_SOL]->GetNodes()->GetVelocity2(iPoint);

      /*--- Compute projections ---*/
      ProjVel = 0.0; bcn = 0.0; vn = 0.0; phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel -= Velocity[iDim]*Normal[iDim];
        bcn     += d[iDim]*UnitNormal[iDim];
        vn      += Velocity[iDim]*UnitNormal[iDim];
        phin    += Psi[iDim+1]*UnitNormal[iDim];
      }

      /*--- Extra boundary term for grid movement ---*/
      if (grid_movement) {
        su2double ProjGridVel = 0.0;
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        phin -= Psi[nVar-1]*ProjGridVel;
      }

      /*--- Introduce the boundary condition ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Psi[iDim+1] -= ( phin - bcn ) * UnitNormal[iDim];

      /*--- Inner products after introducing BC (Psi has changed) ---*/
      phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
      for (iDim = 0; iDim < nDim; iDim++) {
        phis1 -= Normal[iDim]*Psi[iDim+1];
        phis2 += Velocity[iDim]*Psi[iDim+1];
      }

      /*--- Flux of the Euler wall ---*/
      Residual[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual[iDim+1] = ProjVel * Psi[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
      Residual[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One;

      /*--- Flux adjustment for grid movement ---*/
      if (grid_movement) {
        su2double ProjGridVel = 0.0;
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel -= GridVel[iDim]*Normal[iDim];
        Residual[0] -= ProjGridVel*Psi[0];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] -= ProjGridVel*Psi[iDim+1];
        Residual[nVar-1] -= ProjGridVel*Psi[nVar-1];
      }

      if (implicit) {

        /*--- Adjoint density ---*/
        Jacobian_ii[0][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[0][iDim+1] = -ProjVel * (Velocity[iDim] - UnitNormal[iDim] * vn);
        Jacobian_ii[0][nVar-1] = -ProjVel * Enthalpy;

        /*--- Adjoint velocities ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          Jacobian_ii[iDim+1][0] = -Normal[iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_ii[iDim+1][jDim+1] = -ProjVel*(UnitNormal[jDim]*UnitNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitNormal[jDim] * vn));
          Jacobian_ii[iDim+1][iDim+1] += ProjVel;
          Jacobian_ii[iDim+1][nVar-1] = -Normal[iDim] * Enthalpy;
        }

        /*--- Adjoint energy ---*/
        Jacobian_ii[nVar-1][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[nVar-1][iDim+1] = 0.0;
        Jacobian_ii[nVar-1][nVar-1] = ProjVel;

        /*--- Jacobian contribution due to grid movement ---*/
        if (grid_movement) {
          su2double ProjGridVel = 0.0;
          su2double *GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel -= GridVel[iDim]*Normal[iDim];
          Jacobian_ii[0][0] -= ProjGridVel;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[iDim+1][iDim+1] -= ProjGridVel;
          Jacobian_ii[nVar-1][nVar-1] -= ProjGridVel;
        }

        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }

      /*--- Update residual ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      
    }
  }
  
  delete [] Velocity;
  delete [] UnitNormal;
  delete [] Psi;
  
}

void CAdjEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint;
  su2double *Normal, ProjVel = 0.0, vn = 0.0, Area, *UnitNormal,
  *Psi_domain, *Psi_sym, *Velocity, Enthalpy = 0.0,
  sq_vel, phin, phis1, phis2;
  unsigned short iDim, iVar, jDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();

  Normal = new su2double[nDim];
  UnitNormal = new su2double[nDim];
  Velocity = new su2double[nDim];
  Psi_domain = new su2double[nVar];
  Psi_sym = new su2double[nVar];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      
      Area = 0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim]   = -Normal[iDim]/Area;

      /*--- Create a copy of the adjoint solution ---*/

      for (iVar = 0; iVar < nVar; iVar++) Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);

      /*--- Retrieve flow variables ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim);

      Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);
      sq_vel   = 0.5*solver_container[FLOW_SOL]->GetNodes()->GetVelocity2(iPoint);

      /*--- Compute projections ---*/

      ProjVel = 0.0; vn = 0.0; phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel -= Velocity[iDim]*Normal[iDim];
        vn      += Velocity[iDim]*UnitNormal[iDim];
        phin    += Psi_domain[iDim+1]*UnitNormal[iDim];
      }

      /*--- Grid Movement ---*/

      if (grid_movement) {
        su2double ProjGridVel = 0.0;
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        }
        phin -= Psi_domain[nVar-1]*ProjGridVel;
      }

      /*--- Introduce the boundary condition ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Psi_domain[iDim+1] -= phin * UnitNormal[iDim];

      /*--- Inner products after introducing BC (Psi has changed) ---*/

      phis1 = 0.0; phis2 = Psi_domain[0] + Enthalpy * Psi_domain[nVar-1];
      for (iDim = 0; iDim < nDim; iDim++) {
        phis1 -= Normal[iDim]*Psi_domain[iDim+1];
        phis2 += Velocity[iDim]*Psi_domain[iDim+1];
      }

      /*--- Flux of the Euler wall ---*/

      Residual_i[0] = ProjVel * Psi_domain[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_i[iDim+1] = ProjVel * Psi_domain[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
      Residual_i[nVar-1] = ProjVel * Psi_domain[nVar-1] + phis1 * Gamma_Minus_One;

      /*--- Grid Movement ---*/

      if (grid_movement) {
        su2double ProjGridVel = 0.0;
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel -= GridVel[iDim]*Normal[iDim];
        Residual_i[0] -= ProjGridVel*Psi_domain[0];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_i[iDim+1] -= ProjGridVel*Psi_domain[iDim+1];
        Residual_i[nVar-1] -= ProjGridVel*Psi_domain[nVar-1];
      }
      
      /*--- Update residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit stuff ---*/
      
      if (implicit) {
        
        /*--- Adjoint density ---*/

        Jacobian_ii[0][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[0][iDim+1] = -ProjVel * (Velocity[iDim] - UnitNormal[iDim] * vn);
        Jacobian_ii[0][nVar-1] = -ProjVel * Enthalpy;

        /*--- Adjoint velocities ---*/

        for (iDim = 0; iDim < nDim; iDim++) {
          Jacobian_ii[iDim+1][0] = -Normal[iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_ii[iDim+1][jDim+1] = -ProjVel*(UnitNormal[jDim]*UnitNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitNormal[jDim] * vn));
          Jacobian_ii[iDim+1][iDim+1] += ProjVel;
          Jacobian_ii[iDim+1][nVar-1] = -Normal[iDim] * Enthalpy;
        }

        /*--- Adjoint energy ---*/

        Jacobian_ii[nVar-1][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[nVar-1][iDim+1] = 0.0;
        Jacobian_ii[nVar-1][nVar-1] = ProjVel;

        /*--- Contribution from grid movement ---*/

        if (grid_movement) {
          su2double ProjGridVel = 0.0;
          su2double *GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel -= GridVel[iDim]*Normal[iDim];
          Jacobian_ii[0][0] -= ProjGridVel;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[iDim+1][iDim+1] -= ProjGridVel;
          Jacobian_ii[nVar-1][nVar-1] -= ProjGridVel;
        }
        
         /*--- Update jacobian ---*/
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        
      }

    }
  }
  
  delete [] Velocity;
  delete [] Psi_domain;
  delete [] Psi_sym;
  delete [] Normal;
  delete [] UnitNormal;
  
}

void CAdjEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                            CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
  unsigned short iDim, iVar;
  su2double *V_i;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *Psi_i = new su2double[nVar];
  su2double *Psi_j = new su2double[nVar];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
      
      /*--- Store the solution for both points ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Psi_j[iVar] = GetDonorAdjVar(val_marker, iVertex, iVar);
      }
      
      /*--- Set adjoint Variables ---*/
      
      numerics->SetAdjointVar(Psi_i, Psi_j);
      
      /*--- Conservative variables w/o reconstruction (the same at both points) ---*/
      
      V_i = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      numerics->SetPrimitive(V_i, V_i);
      
      /*--- Set Normal ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      numerics->SetNormal(Normal);
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add Residuals and Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
      if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
  }
  
  
  delete[] Normal;
  delete[] Psi_i;
  delete[] Psi_j;
  
}

void CAdjEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                            CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
  unsigned short iDim, iVar;
  su2double *V_i, *IntBoundary_Jump;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *Psi_i = new su2double[nVar];
  su2double *Psi_j = new su2double[nVar];
  su2double *Psi_out = new su2double[nVar];
  su2double *Psi_in = new su2double[nVar];
  su2double *Psi_out_ghost = new su2double[nVar];
  su2double *Psi_in_ghost = new su2double[nVar];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
      
      /*--- Store the solution for both points ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Psi_j[iVar] = GetDonorAdjVar(val_marker, iVertex, iVar);
      }
      
      /*--- Set Normal ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      numerics->SetNormal(Normal);
      
      /*--- If equivalent area or nearfield pressure condition ---*/
      
      if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) ||
          (config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {
        
        /*--- Read the jump ---*/
        
        IntBoundary_Jump = nodes->GetIntBoundary_Jump(iPoint);
        
        /*--- Inner point ---*/
        
        if (Normal[nDim-1] < 0.0)  {
          for (iVar = 0; iVar < nVar; iVar++) {
            Psi_in[iVar] = Psi_i[iVar]; Psi_out[iVar] = Psi_j[iVar];
            Psi_in_ghost[iVar] = Psi_out[iVar] - IntBoundary_Jump[iVar];
          }
          numerics->SetAdjointVar(Psi_in, Psi_in_ghost);
        }
        
        /*--- Outer point ---*/
        
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Psi_in[iVar] = Psi_j[iVar]; Psi_out[iVar] = Psi_i[iVar];
            Psi_out_ghost[iVar] =  Psi_in[iVar] + IntBoundary_Jump[iVar];
          }
          numerics->SetAdjointVar(Psi_out, Psi_out_ghost);
        }
      }
      else {
        
        /*--- Just do a periodic BC ---*/
        
        numerics->SetAdjointVar(Psi_i, Psi_j);
        
      }
      
      /*--- Conservative variables w/o reconstruction (the same at both points) ---*/
      
      V_i = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      numerics->SetPrimitive(V_i, V_i);
      
      /*--- Compute residual ---*/
      
      numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add Residuals and Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
      if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
  }
  
  delete[] Normal;
  delete[] Psi_i;
  delete[] Psi_j;
  delete[] Psi_out;
  delete[] Psi_in;
  delete[] Psi_out_ghost;
  delete[] Psi_in_ghost;
  
}

void CAdjEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, Point_Normal;
  unsigned short iVar, iDim;
  su2double *Normal, *V_domain, *V_infty, *Psi_domain, *Psi_infty;
  
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_infty = new su2double[nVar];
  
  /*--- Loop over all the vertices ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- If the node belongs to the domain ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Set the normal vector ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Adjoint flow solution at the wall ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
        Psi_infty[iVar] = 0.0;
      }
      conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the upwind flux ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
      /*--- Viscous residual contribution, it doesn't work ---*/
      
      if (config->GetViscous()) {
        
        /*--- Points in edge, coordinates and normal vector---*/
        
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Normal);
        
        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/
        
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_infty);
        
        /*--- Gradient and limiter of Adjoint Variables ---*/
        
        visc_numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));
        
        /*--- Compute residual ---*/
        
        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
        
        /*--- Update adjoint viscous residual ---*/
        
        LinSysRes.SubtractBlock(iPoint, Residual_i);
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        
      }
      
    }
  }
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_infty;
}

void CAdjEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
    CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain, *Normal, *Psi_domain, *Psi_inlet;

  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Adjoint flow solution at the boundary ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);

      /*--- Construct the flow & adjoint states at the inlet ---*/
      /*--- Supersonic Inlet: All characteristic are exiting: using nearest neighbor to set value ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Psi_inlet[iVar] = Psi_domain[iVar];

      /*--- Set the flow and adjoint states in the solver ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_inlet);

      /*--- Grid Movement ---*/

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {
        
        /*--- Index of the closest interior node ---*/

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Points in edge, coordinates and normal vector---*/

        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Normal);

        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_inlet);

        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_inlet;
  
}

void CAdjEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                          CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain, *Normal, *Psi_domain, *Psi_outlet;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_outlet = new su2double[nVar];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Allocate the value at the inlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
      Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
      
      /*--- Construct the flow & adjoint states at the inlet ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
      Psi_outlet[iVar] = 0.0;
      
      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_outlet);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
      /*--- Viscous residual contribution (check again, Point_Normal was not being initialized before) ---*/

      if (config->GetViscous()) {

        /*--- Index of the closest interior node ---*/
        
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
        
        /*--- Points in edge, coordinates and normal vector---*/

        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Normal);

        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_outlet;
  
}

void CAdjEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Velocity[3], bcn, phin, Area, UnitNormal[3],
  ProjGridVel, *GridVel;
  su2double *V_inlet, *V_domain, *Normal, *Psi_domain, *Psi_inlet;

  unsigned short Kind_Inlet = config->GetKind_Inlet();
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
      
      /*--- Construct the flow & adjoint states at the inlet ---*/

      /*--- Subsonic, compressible inflow: first build the flow state
         using the same method as the direct problem. Then, based on
         those conservative values, compute the characteristic-based
         adjoint boundary condition. The boundary update to be applied
         depends on whether total conditions or mass flow are specified. ---*/

      switch (Kind_Inlet) {

      /*--- Total properties have been specified at the inlet. ---*/
      case TOTAL_CONDITIONS:

        /*--- Adjoint solution at the inlet. Set to zero for now
             but should be replaced with derived expression for this type of
             inlet. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Psi_inlet[iVar] = 0.0;

        break;

        /*--- Mass flow has been specified at the inlet. ---*/
      case MASS_FLOW:

        /*--- Get primitives from current inlet state. ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim);

        /*--- Retrieve current adjoint solution values at the boundary ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Psi_inlet[iVar] = nodes->GetSolution(iPoint,iVar);

        /*--- Some terms needed for the adjoint BC ---*/
        bcn = 0.0; phin = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          bcn  -= (Gamma/Gamma_Minus_One)*Velocity[iDim]*UnitNormal[iDim];
          phin += Psi_domain[iDim+1]*UnitNormal[iDim];
        }

        /*--- Extra boundary term for grid movement ---*/
        if (grid_movement) {
          ProjGridVel = 0.0;
          GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
          bcn -= (1.0/Gamma_Minus_One)*ProjGridVel;
        }

        /*--- Impose value for PsiE based on hand-derived expression. ---*/
        Psi_inlet[nVar-1] = -phin*(1.0/bcn);

        break;
      }

      
      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_inlet);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {
        /*--- Index of the closest interior node ---*/

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Points in edge, coordinates and normal vector---*/

        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Normal);

        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_inlet);

        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_inlet;
  
}

void CAdjEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Pressure=0.0, P_Exit=0.0,  Velocity2 = 0.0, Area=0.0, Density=0.0,
      Vn = 0.0, SoundSpeed = 0.0, Vn_Exit=0.0, ProjGridVel = 0.0,
      Riemann=0.0, Entropy=0.0, Vn_rel=0.0;
  su2double Velocity[3], UnitNormal[3];
  su2double *V_outlet, *V_domain, *Psi_domain, *Psi_outlet, *Normal;
  su2double a1=0.0, a2=0.0; /*Placeholder terms to simplify expressions/ repeated terms*/
  /*Gradient terms for the generalized boundary */
  su2double density_gradient=0.0, pressure_gradient=0.0, velocity_gradient=0.0;

  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  su2double Weight_ObjFunc = 1.0;  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  string Monitoring_Tag;
  unsigned short jMarker=0, iMarker_Monitoring=0;

  Psi_domain = new su2double [nVar]; Psi_outlet = new su2double [nVar];
  Normal = new su2double[nDim];

  /*--- Identify marker monitoring index ---*/
  for (jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
    Monitoring_Tag = config->GetMarker_Monitoring_TagBound(jMarker);
    if (Monitoring_Tag==Marker_Tag)
      iMarker_Monitoring = jMarker;
  }

  Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);

  /*--- Loop over all the vertices ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- If the node belong to the domain ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Set the normal point ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

      /*--- Allocate the value at the outlet ---*/

      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Adjoint flow solution at the boundary ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);

      /*--- Construct the flow & adjoint states at the outlet ---*/

      /*--- Retrieve the specified back pressure for this outlet, Non-dim. the inputs if necessary. ---*/

      P_Exit = config->GetOutlet_Pressure(Marker_Tag);
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/

      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      /*--- Extra boundary term for grid movement ---*/
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
      }

      Pressure = V_domain[nDim+1];
      SoundSpeed = sqrt(Pressure*Gamma/Density);

      /*--- Set Adjoint variables to 0 initially ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_outlet[iVar] = 0.0;
      }

      if (Vn > SoundSpeed) {
        /*--- Objective-dependent additions to energy term ---*/
        Vn_Exit = Vn; /* Vn_Exit comes from Reiman conditions in subsonic case*/
        Vn_rel = Vn_Exit-ProjGridVel;
        /* Repeated term */
        a1 = Gamma_Minus_One/(Vn_rel*Vn_rel-SoundSpeed*SoundSpeed);

        switch (config->GetKind_ObjFunc(iMarker_Monitoring)) {
        case SURFACE_TOTAL_PRESSURE:
          /*--- Total Pressure term. NOTE: this is AREA averaged
           * Additional terms are added later (as they are common between subsonic,
           * supersonic equations) ---*/
          Velocity2  = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          a2 = Pressure*(Gamma/Gamma_Minus_One)*pow((1.0+Gamma_Minus_One*Density*Velocity2/(2.0*Gamma*Pressure)),1.0/Gamma_Minus_One);
          density_gradient = a2*(Gamma_Minus_One*Velocity2/(2.0*Gamma*Pressure));
          velocity_gradient = 0.0;
          for (iDim=0; iDim<nDim; iDim++)
            velocity_gradient+=a2*Gamma_Minus_One*Density/(Gamma*Pressure)*Velocity[iDim]*UnitNormal[iDim];
          pressure_gradient = a2*(-Gamma_Minus_One*Density*Velocity2/(2.0*Gamma*pow(Pressure,2.0)))+pow((1.0+Gamma_Minus_One*Density*Velocity2/(2.0*Gamma*Pressure)),(Gamma/Gamma_Minus_One));
          Psi_outlet[nDim+1]+=Weight_ObjFunc*a1*(density_gradient/Vn_rel+pressure_gradient*Vn_rel-velocity_gradient/Density);
          break;
        case SURFACE_STATIC_PRESSURE:
          /*Area averaged static pressure*/
          /*--- Note: further terms are NOT added later for this case, only energy term is modified ---*/
          Psi_outlet[nDim+1]+=Weight_ObjFunc*(a1*Vn_Exit);
          break;
        default:
          break;
        }
      } else {
        /*---Subsonic Case: Psi-rho E term from volume, objective-specific terms which are common
           * between subsonic and supersonic cases are added later  ---*/
        /*--- Compute Riemann constant ---*/
        Entropy = Pressure*pow(1.0/Density, Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        /*--- Update velocity terms ---*/
        Vn_rel  = Vn_Exit-ProjGridVel;

        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }

        /*--- Extra boundary term for grid movement ---*/

        if (grid_movement) {
          ProjGridVel = 0.0;
          su2double *GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        }

        /*--- Impose values for PsiRho & Phi using PsiE from domain. ---*/

        Psi_outlet[nVar-1] = Psi_domain[nVar-1];

      }

        /*--- When Psi_outlet[nVar-1] is not 0, the other terms of Psi_outlet must be updated
    This occurs when subsonic, or for certain objective functions ---*/
        if ( Psi_outlet[nVar-1] !=0.0 ) {
          /*--- Shorthand for repeated term in the boundary conditions ---*/
          a1 = 0.0;
          if (Vn!=0.0)
            a1 = SoundSpeed*SoundSpeed/Gamma_Minus_One/Vn;
          Psi_outlet[0] += Psi_outlet[nVar-1]*(Velocity2*0.5+Vn_rel*a1);
          for (iDim = 0; iDim < nDim; iDim++) {
            Psi_outlet[iDim+1] += -Psi_outlet[nVar-1]*(a1*UnitNormal[iDim] + Velocity[iDim]);
          }
        }


      /*--- Add terms for objective functions where additions are needed outside the energy term
       *     Terms which are added to the energy term are taken care of in the supersonic section above ---*/
      switch (config->GetKind_ObjFunc(iMarker_Monitoring)) {
      case SURFACE_MASSFLOW:
        Psi_outlet[0]+=Weight_ObjFunc;
        break;
      case SURFACE_TOTAL_PRESSURE:
        /*--- For total pressure objective function. NOTE: this is AREA averaged term---*/
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        if (Vn_Exit !=0.0) {
          a2 = Pressure*(Gamma/Gamma_Minus_One)*pow((1.0+Gamma_Minus_One*Density*Velocity2/(2.0*Gamma*Pressure)),1.0/(Gamma_Minus_One));
          density_gradient = a2*(Gamma_Minus_One*Velocity2/(2.0*Gamma*Pressure));
          velocity_gradient=a2*Gamma_Minus_One*Density/(Gamma*Pressure); // re-using variable as the constant multiplying V[i] for dj/dvi
          Psi_outlet[0]+=Weight_ObjFunc*(density_gradient*2.0/Vn_Exit);
          for (iDim=0; iDim<nDim; iDim++) {
            Psi_outlet[0]-=Weight_ObjFunc*(velocity_gradient*Velocity[iDim]*Velocity[iDim]/(Density*Vn_Exit));
            Psi_outlet[iDim+1] += Weight_ObjFunc*(velocity_gradient*Velocity[iDim]/(Density*Vn_Exit) - UnitNormal[iDim]*density_gradient/(Vn_Exit*Vn_Exit));
          }
        }
        break;
      case SURFACE_STATIC_PRESSURE:
        /*0.0*/
        break;
      default:
        break;
      }

      /*--- Set the flow and adjoint states in the solver ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

      /*--- Grid Movement ---*/

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
          Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {

        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_outlet[nDim+5] = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        V_outlet[nDim+6] = solver_container[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint);

        /*--- Points in edge, coordinates and normal vector---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

        /*--- Turbulent kinetic energy ---*/

        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0), solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_outlet;
  
}

void CAdjEulerSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  su2double *Normal, *V_domain, *V_inflow, *Psi_domain, *Psi_inflow, P_Fan, Velocity[3],
  Velocity2, Density, Vn, UnitNormal[3], Area, a1;
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inflow = new su2double[nVar];
  
  /*--- Loop over all the vertices ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- If the node belong to the domain ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Allocate the value at the inflow ---*/
      
      V_inflow = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
      
      /*--- Subsonic flow is assumed, note that there is no non-dimensionalization. ---*/
      
      P_Fan = config->GetInflow_Pressure(Marker_Tag);
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      
      /*---Subsonic Case(s) using Riemann invariants ---*/
      
      //      Pressure = V_domain[nDim+1];
      //      SoundSpeed = sqrt(Pressure*Gamma/Density);
      //      Mach_Fan  = sqrt(Velocity2)/SoundSpeed;
      //      Entropy = Pressure*pow(1.0/Density, Gamma);
      //      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
      
      //      /*--- Compute (Vn - Ubn).n term for use in the BC.
      //      Compute the new fictious state at the outlet ---*/
      
      //      Density    = pow(P_Fan/Entropy,1.0/Gamma);
      //      SoundSpeed = sqrt(Gamma*P_Fan/Density);
      //      Vn_Fan    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      //      Velocity2  = 0.0;
      //      for (iDim = 0; iDim < nDim; iDim++) {
      //        Velocity[iDim] = Velocity[iDim] + (Vn_Fan-Vn)*UnitNormal[iDim];
      //        Velocity2 += Velocity[iDim]*Velocity[iDim];
      //      }
      
      /*--- Shorthand for repeated term in the boundary conditions ---*/
      
      a1 = sqrt(Gamma*P_Fan/Density)/(Gamma_Minus_One);
      
      /*--- Impose values for PsiRho & Phi using PsiE from domain. ---*/
      
      Psi_inflow[nVar-1] = Psi_domain[nVar-1];
      Psi_inflow[0] = 0.5*Psi_inflow[nVar-1]*Velocity2;
      for (iDim = 0; iDim < nDim; iDim++) {
        Psi_inflow[0]   += Psi_inflow[nVar-1]*a1*Velocity[iDim]*UnitNormal[iDim];
        Psi_inflow[iDim+1] = -Psi_inflow[nVar-1]*(a1*UnitNormal[iDim] + Velocity[iDim]);
      }
      
      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inflow);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_inflow);
      
      /*--- Compute the residual ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain;
  delete [] Psi_inflow;
  
}

void CAdjEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Area, *Normal, *V_domain, *V_exhaust, *Psi_domain, *Psi_exhaust;
  unsigned short iVar, iDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_exhaust = new su2double[nVar];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Allocate the value at the exhaust ---*/
      
      V_exhaust = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
      
      /*--- Adjoint flow solution at the exhaust (this should be improved using characteristics bc) ---*/
      
      Psi_exhaust[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Psi_exhaust[iDim+1] = nodes->GetSolution(Point_Normal,iDim+1);
      Psi_exhaust[nDim+1] = 0.0;
      
      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_exhaust);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_exhaust);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
  }
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_exhaust;
  
}

void CAdjEulerSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  su2double *Normal, *V_domain, *V_inlet, *Psi_domain, *Psi_inlet;
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, GlobalIndex_inlet, GlobalIndex;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];
  
  /*--- Loop over all the vertices ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    GlobalIndex_inlet = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
    
    /*--- If the node belong to the domain ---*/
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex != GlobalIndex_inlet)) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
        Psi_inlet[iVar] =  0.0;// nodes->GetSolution(iPoint,iVar);
      }
      
#ifdef CHECK
      su2double UnitNormal[3], Area=0.0;
      
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Subsonic flow is assumed, note that there is no non-dimensionalization. ---*/
      
      P_Fan = V_domain[nDim+1] ;
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      
      /*---Subsonic Case(s) using Riemann invariants ---*/
      
      Pressure = V_domain[nDim+1];
      SoundSpeed = sqrt(Pressure*Gamma/Density);
      
      Entropy = Pressure*pow(1.0/Density, Gamma);
      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
      
      /*--- Compute (Vn - Ubn).n term for use in the BC.
       Compute the new fictious state at the outlet ---*/
      
      Density    = pow(P_Fan/Entropy,1.0/Gamma);
      SoundSpeed = sqrt(Gamma*P_Fan/Density);
      Vn_Fan    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      
      Velocity2  = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = Velocity[iDim] + (Vn_Fan-Vn)*UnitNormal[iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      
      /*--- Impose values for PsiRho & Phi using PsiE from domain. ---*/
      
      Psi_inlet[nVar-1] = Psi_domain[nVar-1];
      
      a1 = SoundSpeed*SoundSpeed/Gamma_Minus_One/Vn;
      
      Psi_inlet[0] = Psi_inlet[nVar-1]*(Velocity2*0.5+Vn_Fan*a1);
      Psi_inlet[iDim+1] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Psi_inlet[iDim+1] += -Psi_inlet[nVar-1]*(a1*UnitNormal[iDim] + Velocity[iDim]);
      }
      
#endif
      
      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_inlet);
      
      /*--- Compute the residual ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain;
  delete [] Psi_inlet;
  
}

void CAdjEulerSolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_inlet, GlobalIndex;
  su2double *Normal, *V_domain, *V_outlet, *Psi_domain, *Psi_outlet;
  unsigned short iVar, iDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar];
  Psi_outlet = new su2double[nVar];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    GlobalIndex_inlet = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
    
    /*--- Check that the node belongs to the domain (i.e., not a halo node) and to discard the perimeter ---*/
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex != GlobalIndex_inlet)) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = nodes->GetSolution(iPoint,iVar);
        Psi_outlet[iVar] = 0.0; //nodes->GetSolution(iPoint,iVar);
      }
      
#ifdef CHECK
      unsigned long Point_Normal;
      /*--- Index of the closest interior node ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Adjoint flow solution at the outlet (this should be improved using characteristics bc) ---*/
      
      Psi_outlet[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Psi_outlet[iDim+1] = nodes->GetSolution(Point_Normal,iDim+1);
      Psi_outlet[nDim+1] = 0.0;
      
#endif
      
      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_outlet);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
    
  }
  
  delete [] Normal;
  delete [] Psi_domain;
  delete [] Psi_outlet;
  
}

void CAdjEulerSolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker, bool val_inlet_surface) {
  
}

void CAdjEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                           unsigned short iMesh, unsigned short RunTime_EqSystem) {
  unsigned short iVar, jVar;
  unsigned long iPoint;
  su2double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool Grid_Movement = config->GetGrid_Movement();
  
  /*--- loop over points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Solution at time n-1, n and n+1 ---*/
    U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
    U_time_n   = nodes->GetSolution_time_n(iPoint);
    U_time_nP1 = nodes->GetSolution(iPoint);
    
    /*--- Volume at time n-1 and n ---*/
    if (Grid_Movement) {
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_n = geometry->node[iPoint]->GetVolume_n();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
    }
    else {
      Volume_nM1 = geometry->node[iPoint]->GetVolume();
      Volume_n = geometry->node[iPoint]->GetVolume();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
    }
    
    /*--- Time Step ---*/
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute Residual ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
      if (config->GetTime_Marching() == DT_STEPPING_2ND)
        Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
                          +  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
    }
        
    /*--- Add Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual);
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
        
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
      }
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
  
}

void CAdjEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  string filename, restart_filename;

  /*--- Restart the solution from file information ---*/

  filename         = config->GetSolution_AdjFileName();
  restart_filename = config->GetObjFunc_Extension(filename);
  restart_filename = config->GetFilename(restart_filename, "", val_iter);

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];

      nodes->SetSolution(iPoint_Local,Solution);
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We also call the preprocessing routine
   on the fine level in order to have all necessary quantities updated. ---*/
  
  solver[MESH_0][ADJFLOW_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][ADJFLOW_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][ADJFLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][ADJFLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][ADJFLOW_SOL]->GetNodes()->SetSolution(iPoint,Solution);
    }
    solver[iMesh][ADJFLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][ADJFLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][ADJFLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

CAdjNSSolver::CAdjNSSolver(void) : CAdjEulerSolver() { }

CAdjNSSolver::CAdjNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CAdjEulerSolver() {
  unsigned long iPoint, iVertex;
  string text_line, mesh_filename;
  unsigned short iDim, iVar, iMarker, nLineLets;
  ifstream restart_file;
  string filename, AdjExt;

  su2double RefArea    = config->GetRefArea();
  su2double RefDensity  = config->GetDensity_FreeStreamND();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double Mach_Motion     = config->GetMach_Motion();
  su2double Area=0.0, *Normal = NULL, myArea_Monitored;
  su2double RefVel2, Mach2Vel, Weight_ObjFunc, factor;
  su2double *Velocity_Inf;

  string Marker_Tag, Monitoring_Tag;
  unsigned short iMarker_Monitoring, jMarker, ObjFunc;
  bool grid_movement  = config->GetGrid_Movement();
  bool restart = config->GetRestart();
  
  /*--- Norm heat flux objective test ---*/
  pnorm = 1.0;
  if (config->GetKind_ObjFunc()==MAXIMUM_HEATFLUX)
    pnorm = 8.0; // Matches MaxNorm defined in solver_direct_mean.
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim         = geometry->GetnDim();
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  nVar = nDim + 2;

  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define some auxiliary arrays related to the residual ---*/
  
  Point_Max    = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]     = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]   = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]   = 0.0;
  Res_Conv_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]   = 0.0;
  Res_Visc_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]   = 0.0;
  Res_Conv_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]   = 0.0;
  Res_Visc_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary arrays related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

  /*--- Define some auxiliary arrays related to the flow solution ---*/
  
  FlowPrimVar_i = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_i[iVar] = 0.0;
  FlowPrimVar_j = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Solution and residual vectors ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
    Jacobian_ii = new su2double*[nVar];
    Jacobian_ij = new su2double*[nVar];
    Jacobian_ji = new su2double*[nVar];
    Jacobian_jj = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_ii[iVar] = new su2double[nVar];
      Jacobian_ij[iVar] = new su2double[nVar];
      Jacobian_ji[iVar] = new su2double[nVar];
      Jacobian_jj[iVar] = new su2double[nVar];
    }
    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Array structures for computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Sensitivity definition and coefficient on all markers ---*/
  CSensitivity = new su2double* [nMarker];
  for (iMarker=0; iMarker<nMarker; iMarker++) {
    CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]];
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  DonorAdjVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorAdjVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorAdjVar[iMarker][iVertex] = new su2double [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        DonorAdjVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }
  
  Sens_Geo   = new su2double[nMarker];
  Sens_Mach  = new su2double[nMarker];
  Sens_AoA   = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];
  Sens_BPress = new su2double[nMarker];
  
  /*--- Initialize sensitivities to zero ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker]   = 0.0;
    Sens_Mach[iMarker]  = 0.0;
    Sens_AoA[iMarker]   = 0.0;
    Sens_Press[iMarker] = 0.0;
    Sens_Temp[iMarker]  = 0.0;
    Sens_BPress[iMarker] = 0.0;
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
      CSensitivity[iMarker][iVertex] = 0.0;
  }
  
  /*--- Initialize the adjoint variables to zero (infinity state) ---*/
  PsiRho_Inf = 0.0;
  if ((config->GetKind_ObjFunc() == TOTAL_HEATFLUX) ||
      (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX) ||
      (config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX))
    PsiE_Inf = 1.0;
  else
    PsiE_Inf = 0.0;
  Phi_Inf = new su2double [nDim];
  Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
  if (nDim == 3) Phi_Inf[2] = 0.0;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CAdjNSVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Read the restart metadata. ---*/

  if (restart && (iMesh == MESH_0)) {
    mesh_filename = config->GetSolution_AdjFileName();
    filename      = config->GetObjFunc_Extension(mesh_filename);
//    Read_SU2_Restart_Metadata(geometry, config, true, filename);
  }

  /*--- Calculate area monitored for area-averaged-outflow-quantity-based objectives ---*/
   myArea_Monitored = 0.0;
   for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
     if (config->GetKind_ObjFunc(iMarker_Monitoring)==SURFACE_TOTAL_PRESSURE ||
         config->GetKind_ObjFunc(iMarker_Monitoring)==SURFACE_STATIC_PRESSURE) {

       Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
       /*-- Find the marker index ---*/
       iMarker = 0;
       for (jMarker= 0; jMarker < config->GetnMarker_All(); jMarker++) {
         Marker_Tag = config->GetMarker_All_TagBound(jMarker);
         if (Marker_Tag == Monitoring_Tag) {
           iMarker = jMarker;
           for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
             iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
             if (geometry->node[iPoint]->GetDomain()) {
               Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
               Area = 0.0;
               for (iDim = 0; iDim < nDim; iDim++)
                 Area += Normal[iDim]*Normal[iDim];
               myArea_Monitored += sqrt (Area);
             }
           }
           break;
         }
       }
     }
   }


 #ifdef HAVE_MPI
   Area_Monitored = 0.0;
   SU2_MPI::Allreduce(&myArea_Monitored, &Area_Monitored, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 #else
   Area_Monitored = myArea_Monitored;
 #endif

   if (config->GetnObj() > 1 && iMesh == MESH_0) {
     if (grid_movement) {
       Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
       RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
     }
     else {
       Velocity_Inf = config->GetVelocity_FreeStreamND();
       RefVel2 = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
     }

     /*--- Objective scaling: a factor must be applied to certain objectives ---*/
     for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
         Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);

         factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);

         ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);
         if ((ObjFunc == INVERSE_DESIGN_HEATFLUX) ||
             (ObjFunc == TOTAL_HEATFLUX) || (ObjFunc == MAXIMUM_HEATFLUX) ||
             (ObjFunc == SURFACE_MASSFLOW) ) factor = 1.0;

        if ((ObjFunc == SURFACE_TOTAL_PRESSURE) || (ObjFunc == SURFACE_STATIC_PRESSURE))
          factor = 1.0/Area_Monitored;

        Weight_ObjFunc = Weight_ObjFunc*factor;
        config->SetWeight_ObjFunc(iMarker_Monitoring, Weight_ObjFunc);

     }
   }

  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
}

CAdjNSSolver::~CAdjNSSolver(void) {
  
}


void CAdjNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh, unsigned long Iteration) {

  /*--- Use the flow solution to update the time step
   *    The time step depends on the characteristic velocity, which is the same
   *    for the adjoint and flow solutions, albeit in the opposite direction. ---*/
  solver_container[FLOW_SOL]->SetTime_Step(geometry, solver_container, config, iMesh, Iteration);
  
}

void CAdjNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double SharpEdge_Distance;
  bool RightSol = true;
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool limiter        = (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);
  bool fixed_cl       = config->GetFixed_CL_Mode();
  bool eval_dof_dcx   = config->GetEval_dOF_dCX();

  /*--- Update the objective function coefficient to guarantee zero gradient. ---*/
  
  if (fixed_cl && eval_dof_dcx) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Residual initialization ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Get the distance form a sharp edge ---*/
    
    SharpEdge_Distance = geometry->node[iPoint]->GetSharpEdge_Distance();
    
    /*--- Initialize the non-physical points vector ---*/
    
    nodes->SetNon_Physical(iPoint,false);
    
    /*--- Set the primitive variables compressible
     adjoint variables ---*/
    
    RightSol = nodes->SetPrimVar(iPoint,SharpEdge_Distance, false, config);
    if (!RightSol) { nodes->SetNon_Physical(iPoint,true); ErrorCounter++; }
    
    /*--- Initialize the convective residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Compute gradients adj for solution reconstruction and viscous term ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Limiter computation (upwind reconstruction) ---*/
  
  if (limiter && !Output) SetSolution_Limiter(geometry, config);

  /*--- Compute gradients adj for viscous term coupling ---*/

  if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solver_container[ADJTURB_SOL]->SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solver_container[ADJTURB_SOL]->SetSolution_Gradient_LS(geometry, config);
  }
  
  /*--- Artificial dissipation for centered schemes ---*/
  
  if (center_jst && (iMesh == MESH_0)) {
    SetCentered_Dissipation_Sensor(geometry, config);
    SetUndivided_Laplacian(geometry, config);
  }
  
  /*--- Initialize the Jacobian for implicit integration ---*/
  
  if (implicit) Jacobian.SetValZero();
  
  /*--- Error message ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CAdjNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, coordinates and normal vector---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables w/o reconstruction and adjoint variables w/o reconstruction---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint),
                           solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(jPoint));
    
    numerics->SetAdjointVar(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    
    /*--- Gradient and limiter of Adjoint Variables ---*/
    
    numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(jPoint));
    
    /*--- Compute residual ---*/
    
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

    /*--- Update adjoint viscous residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
  
}

void CAdjNSSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  
  /*--- Loop over all the points, note that we are supposing that primitive and
   adjoint gradients have been computed previously ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Primitive variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint), NULL);

    /*--- Gradient of adjoint variables ---*/
    
    numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), NULL);

    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/
    
    if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
      
      /*--- Turbulent variables w/o reconstruction and its gradient ---*/
      
      numerics->SetTurbVar(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint), NULL);
      
      numerics->SetTurbVarGradient(solver_container[TURB_SOL]->GetNodes()->GetGradient(iPoint), NULL);
      
      /*--- Turbulent adjoint variables w/o reconstruction and its gradient ---*/
      
      numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->GetNodes()->GetSolution(iPoint), NULL);
      
      numerics->SetTurbAdjointGradient(solver_container[ADJTURB_SOL]->GetNodes()->GetGradient(iPoint), NULL);
      
      /*--- Set distance to the surface ---*/
      
      numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
      
    }
    
    /*--- Compute residual ---*/
    
    numerics->ComputeResidual(Residual, config);
    
    /*--- Add to the residual ---*/
    
    LinSysRes.AddBlock(iPoint, Residual);
    
  }
  
  /*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/
  
  if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Points in edge, and normal vector ---*/

      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      second_numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

      /*--- Conservative variables w/o reconstruction ---*/

      second_numerics->SetConservative(solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint),
                                     solver_container[FLOW_SOL]->GetNodes()->GetSolution(jPoint));

      /*--- Gradient of primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint),
                                        solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(jPoint));

      /*--- Viscosity ---*/

      second_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint),
                                         solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint));

      /*--- Turbulent variables w/o reconstruction ---*/

      second_numerics->SetTurbVar(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint),
                                solver_container[TURB_SOL]->GetNodes()->GetSolution(jPoint));

      /*--- Turbulent adjoint variables w/o reconstruction ---*/

      second_numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->GetNodes()->GetSolution(iPoint),
                                       solver_container[ADJTURB_SOL]->GetNodes()->GetSolution(jPoint));

      /*--- Set distance to the surface ---*/

      second_numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), geometry->node[jPoint]->GetWall_Distance());
      
      /*--- Update adjoint viscous residual ---*/
      
      second_numerics->ComputeResidual(Residual, config);
      
      LinSysRes.AddBlock(iPoint, Residual);
      LinSysRes.SubtractBlock(jPoint, Residual);
    }

  }
  
  // WARNING: The rotating frame source term has been placed in the second
  // source term container since the section below is commented. This needs a
  // permanent fix asap!
  
  if (rotating_frame) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the adjoint variables ---*/
      second_numerics->SetAdjointVar(nodes->GetSolution(iPoint),
                                     nodes->GetSolution(iPoint));
      
      /*--- Load the volume of the dual mesh cell ---*/
      second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the adjoint rotating frame source residual ---*/
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
}

void CAdjNSSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, jDim, iMarker, iPos, jPos;
  su2double *d = NULL, **PsiVar_Grad = NULL, **PrimVar_Grad = NULL, div_phi, *Normal = NULL, Area,
  normal_grad_psi5, normal_grad_T, sigma_partial, Laminar_Viscosity = 0.0, heat_flux_factor, temp_sens = 0.0, *Psi = NULL, *U = NULL, Enthalpy, **GridVel_Grad, gradPsi5_v, psi5_tau_partial, psi5_tau_grad_vel, source_v_1, Density, Pressure = 0.0, div_vel, val_turb_ke, vartheta, vartheta_partial, psi5_p_div_vel, Omega[3], rho_v[3] = {0.0,0.0,0.0}, CrossProduct[3], delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}, r, ru, rv, rw, rE, p, T, dp_dr, dp_dru, dp_drv, dp_drw, dp_drE, dH_dr, dH_dru, dH_drv, dH_drw, dH_drE, H, D[3][3], Dd[3], Mach_Inf, eps, scale = 1.0;
  su2double RefVel2, RefDensity, Mach2Vel, *Velocity_Inf, factor;

  su2double *USens = new su2double[nVar];
  su2double *UnitNormal = new su2double[nDim];
  su2double *normal_grad_vel = new su2double[nDim];
  su2double *tang_deriv_psi5 = new su2double[nDim];
  su2double *tang_deriv_T = new su2double[nDim];
  su2double **Sigma = new su2double* [nDim];
  
  for (iDim = 0; iDim < nDim; iDim++)
    Sigma[iDim] = new su2double [nDim];
  
  su2double *normal_grad_gridvel = new su2double[nDim];
  su2double *normal_grad_v_ux =new su2double[nDim];
  su2double **Sigma_Psi5v = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Sigma_Psi5v[iDim] = new su2double [nDim];
  su2double **tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    tau[iDim] = new su2double [nDim];
  su2double *Velocity = new su2double[nDim];
  
  bool rotating_frame    = config->GetRotating_Frame();
  bool grid_movement     = config->GetGrid_Movement();
  su2double RefArea    = config->GetRefArea();
  su2double Mach_Motion     = config->GetMach_Motion();
  unsigned short ObjFunc = config->GetKind_ObjFunc();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double Cp              = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Prandtl_Lam     = config->GetPrandtl_Lam();
  
  if (config->GetSystemMeasurements() == US) scale = 1.0/12.0;
  else scale = 1.0;
  
  /*--- Compute non-dimensional factor. For dynamic meshes, use the motion Mach
   number as a reference value for computing the force coefficients.
   Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  RefDensity  = config->GetDensity_FreeStreamND();
  
  factor = 1.0;
  /*-- For multi-objective problems these scaling factors are applied before solution ---*/
  if (config->GetnObj()==1) {
    factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);

    if ((ObjFunc == INVERSE_DESIGN_HEATFLUX) ||
        (ObjFunc == TOTAL_HEATFLUX) || (ObjFunc == MAXIMUM_HEATFLUX) ||
        (ObjFunc == SURFACE_MASSFLOW))
      factor = 1.0;

    if ((ObjFunc == SURFACE_TOTAL_PRESSURE) || (ObjFunc == SURFACE_STATIC_PRESSURE))
      factor = 1.0/Area_Monitored;

  }


  /*--- Compute gradient of the grid velocity, if applicable ---*/
  
  if (grid_movement)
    SetGridVel_Gradient(geometry, config);
  
  Total_Sens_Geo = 0.0;
  Total_Sens_Mach = 0.0;
  Total_Sens_AoA = 0.0;
  Total_Sens_Press = 0.0;
  Total_Sens_Temp = 0.0;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Sens_Geo[iMarker] = 0.0;
    
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          PsiVar_Grad = nodes->GetGradient(iPoint);
          PrimVar_Grad = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
          
          Laminar_Viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
          
          heat_flux_factor = Cp * Laminar_Viscosity / Prandtl_Lam;
          
          /*--- Compute face area and the unit normal to the surface ---*/
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Normal[iDim]*Normal[iDim]; } Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) { UnitNormal[iDim] = Normal[iDim] / Area; }
          
          /*--- Compute the sensitivity related to the temperature ---*/
          

          normal_grad_psi5 = 0.0; normal_grad_T = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_grad_psi5 += PsiVar_Grad[nVar-1][iDim]*UnitNormal[iDim];
            normal_grad_T += PrimVar_Grad[0][iDim]*UnitNormal[iDim];
          }

          temp_sens = 0.0;
          if (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {

            /*--- Heat Flux Term: temp_sens = (\partial_tg \psi_5)\cdot (k \partial_tg T) ---*/

            for (iDim = 0; iDim < nDim; iDim++) {
              tang_deriv_psi5[iDim] = PsiVar_Grad[nVar-1][iDim] - normal_grad_psi5*UnitNormal[iDim];
              tang_deriv_T[iDim] = PrimVar_Grad[0][iDim] - normal_grad_T*UnitNormal[iDim];
            }
            for (iDim = 0; iDim < nDim; iDim++)
              temp_sens += heat_flux_factor * tang_deriv_psi5[iDim] * tang_deriv_T[iDim];

          } else if (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {

            /*--- Isothermal Term: temp_sens = - k * \partial_n(\psi_5) * \partial_n(T) ---*/

            temp_sens = - heat_flux_factor * normal_grad_psi5 * normal_grad_T;

          }

          
          /*--- Term: sigma_partial = \Sigma_{ji} n_i \partial_n v_j ---*/
          
          div_phi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            div_phi += PsiVar_Grad[iDim+1][iDim];
            for (jDim = 0; jDim < nDim; jDim++)
              Sigma[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
          }
          for (iDim = 0; iDim < nDim; iDim++)
            Sigma[iDim][iDim] -= TWO3*Laminar_Viscosity * div_phi;

          
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_grad_vel[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              normal_grad_vel[iDim] += PrimVar_Grad[iDim+1][jDim]*UnitNormal[jDim];
          }
          
          sigma_partial = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            for (jDim = 0; jDim < nDim; jDim++)
              sigma_partial += UnitNormal[iDim]*Sigma[iDim][jDim]*normal_grad_vel[jDim];
          
          /*--- Compute additional terms in the surface sensitivity for
           moving walls in a rotating frame or dynamic mesh problem. ---*/
          
          if (grid_movement) {
            
            Psi = nodes->GetSolution(iPoint);
            U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
            Density = U[0];
            Pressure = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);
            Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);
            
            /*--- Turbulent kinetic energy ---*/
            
            if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
              val_turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
            else
              val_turb_ke = 0.0;
            
            div_vel = 0.0;
            for (iDim = 0 ; iDim < nDim; iDim++) {
              Velocity[iDim] = U[iDim+1]/Density;
              div_vel += PrimVar_Grad[iDim+1][iDim];
            }
            
            for (iDim = 0 ; iDim < nDim; iDim++)
              for (jDim = 0 ; jDim < nDim; jDim++)
                tau[iDim][jDim] = Laminar_Viscosity*(PrimVar_Grad[jDim+1][iDim] + PrimVar_Grad[iDim+1][jDim])
                - TWO3*Laminar_Viscosity*div_vel*delta[iDim][jDim]
                - TWO3*Density*val_turb_ke*delta[iDim][jDim];
            
            /*--- Form normal_grad_gridvel = \partial_n (u_omega) ---*/
            
            GridVel_Grad = geometry->node[iPoint]->GetGridVel_Grad();
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_grad_gridvel[iDim] = 0.0;
              for (jDim = 0; jDim < nDim; jDim++)
                normal_grad_gridvel[iDim] += GridVel_Grad[iDim][jDim]*UnitNormal[jDim];
            }
            
            /*--- Form normal_grad_v_ux = \partial_n (v - u_omega) ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              normal_grad_v_ux[iDim] = normal_grad_vel[iDim] - normal_grad_gridvel[iDim];
            
            /*--- Form Sigma_Psi5v ---*/
            
            gradPsi5_v = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              gradPsi5_v += PsiVar_Grad[nDim+1][iDim]*Velocity[iDim];
              for (jDim = 0; jDim < nDim; jDim++)
                Sigma_Psi5v[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[nDim+1][iDim]*Velocity[jDim]+PsiVar_Grad[nDim+1][jDim]*Velocity[iDim]);
            }
            for (iDim = 0; iDim < nDim; iDim++)
              Sigma_Psi5v[iDim][iDim] -= TWO3*Laminar_Viscosity * gradPsi5_v;
            
            
            /*--- Now compute terms of the surface sensitivity ---*/
            
            /*--- Form vartheta_partial = \vartheta * \partial_n (v - u_x) . n ---*/
            vartheta = Density*Psi[0] + Density*Enthalpy*Psi[nDim+1];
            for (iDim = 0; iDim < nDim; iDim++) {
              vartheta += U[iDim+1]*Psi[iDim+1];
            }
            vartheta_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              vartheta_partial += vartheta * normal_grad_v_ux[iDim] * UnitNormal[iDim];
            
            /*--- Form sigma_partial = n_i ( \Sigma_Phi_{ij} + \Sigma_Psi5v_{ij} ) \partial_n (v - u_x)_j ---*/
            
            sigma_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                sigma_partial += UnitNormal[iDim]*(Sigma[iDim][jDim]+Sigma_Psi5v[iDim][jDim])*normal_grad_v_ux[jDim];
            
            /*--- Form psi5_tau_partial = \Psi_5 * \partial_n (v - u_x)_i * tau_{ij} * n_j ---*/
            
            psi5_tau_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                psi5_tau_partial -= Psi[nDim+1]*normal_grad_v_ux[iDim]*tau[iDim][jDim]*UnitNormal[jDim];
            
            /*--- Form psi5_p_div_vel = ---*/
            
            psi5_p_div_vel = -Psi[nDim+1]*Pressure*div_vel;
            
            /*--- Form psi5_tau_grad_vel = \Psi_5 * tau_{ij} : \nabla( v ) ---*/
            
            psi5_tau_grad_vel = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                psi5_tau_grad_vel += Psi[nDim+1]*tau[iDim][jDim]*PrimVar_Grad[iDim+1][jDim];
            
            /*--- Retrieve the angular velocity vector ---*/
            
            source_v_1 = 0.0;
            if (rotating_frame) {
              
              for (iDim = 0; iDim < 3; iDim++){
                Omega[iDim] = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
              }
              
              /*--- Calculate momentum source terms as: rho * ( Omega X V ) ---*/
              
              for (iDim = 0; iDim < nDim; iDim++)
                rho_v[iDim] = U[iDim+1];
              if (nDim == 2) rho_v[2] = 0.0;
              
              CrossProduct[0] = Omega[1]*rho_v[2] - Omega[2]*rho_v[1];
              CrossProduct[1] = Omega[2]*rho_v[0] - Omega[0]*rho_v[2];
              CrossProduct[2] = Omega[0]*rho_v[1] - Omega[1]*rho_v[0];
              
              
              for (iDim = 0; iDim < nDim; iDim++) {
                source_v_1 += Psi[iDim+1]*CrossProduct[iDim];
              }
            }
            
            /*--- For simplicity, store all additional terms within sigma_partial ---*/
            
            sigma_partial = sigma_partial + vartheta_partial + psi5_tau_partial + psi5_p_div_vel + psi5_tau_grad_vel + source_v_1;
            
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          
          CSensitivity[iMarker][iVertex] = (sigma_partial - temp_sens) * Area * scale * factor;
            
          /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
          
          if (config->GetSens_Remove_Sharp()) {
            eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
            if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
              CSensitivity[iMarker][iVertex] = 0.0;
          }
          
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex];
          
        }
      }
      
      Total_Sens_Geo += Sens_Geo[iMarker];
      
    }
  }
  
  /*--- Farfield Sensitivity (Mach, AoA, Press, Temp), only for compressible flows ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FAR_FIELD || config->GetMarker_All_KindBC(iMarker) == INLET_FLOW ||
        config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_INLET || config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_OUTLET ||
        config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW  ) {

      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      Sens_BPress[iMarker] = 0.0;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          r = U[0]; ru = U[1]; rv = U[2];
          if (nDim == 2) { rw = 0.0; rE = U[3]; }
          else { rw = U[3]; rE = U[4]; }
          p = Gamma_Minus_One*(rE-(ru*ru + rv*rv + rw*rw)/(2*r));

          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

          H = (rE + p)/r;

          dp_dr = Gamma_Minus_One*(ru*ru + rv*rv + rw*rw)/(2*r*r);
          dp_dru = -Gamma_Minus_One*ru/r;
          dp_drv = -Gamma_Minus_One*rv/r;
          if (nDim == 2) { dp_drw = 0.0; dp_drE = Gamma_Minus_One; }
          else { dp_drw = -Gamma_Minus_One*rw/r; dp_drE = Gamma_Minus_One; }

          dH_dr = (-H + dp_dr)/r; dH_dru = dp_dru/r; dH_drv = dp_drv/r;
          if (nDim == 2) { dH_drw = 0.0; dH_drE = (1 + dp_drE)/r; }
          else { dH_drw = dp_drw/r; dH_drE = (1 + dp_drE)/r; }

          if (nDim == 2) {
            Jacobian_j[0][0] = 0.0;
            Jacobian_j[1][0] = Area*UnitNormal[0];
            Jacobian_j[2][0] = Area*UnitNormal[1];
            Jacobian_j[3][0] = 0.0;

            Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitNormal[0] + (-(ru*rv)/(r*r))*Area*UnitNormal[1];
            Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1];
            Jacobian_j[2][1] = (dp_drv)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[1];
            Jacobian_j[3][1] = (dp_drE)*Area*UnitNormal[0];

            Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitNormal[0] + (-(rv*rv)/(r*r) + dp_dr)*Area*UnitNormal[1];
            Jacobian_j[1][2] = (rv/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[1];
            Jacobian_j[2][2] = (ru/r)*Area*UnitNormal[0] + (2*rv/r + dp_drv)*Area*UnitNormal[1];
            Jacobian_j[3][2] = (dp_drE)*Area*UnitNormal[1];

            Jacobian_j[0][3] = (ru*dH_dr)*Area*UnitNormal[0] + (rv*dH_dr)*Area*UnitNormal[1];
            Jacobian_j[1][3] = (H + ru*dH_dru)*Area*UnitNormal[0] + (rv*dH_dru)*Area*UnitNormal[1];
            Jacobian_j[2][3] = (ru*dH_drv)*Area*UnitNormal[0] + (H + rv*dH_drv)*Area*UnitNormal[1];
            Jacobian_j[3][3] = (ru*dH_drE)*Area*UnitNormal[0] + (rv*dH_drE)*Area*UnitNormal[1];
          }
          else {
            Jacobian_j[0][0] = 0.0;
            Jacobian_j[1][0] = Area*UnitNormal[0];
            Jacobian_j[2][0] = Area*UnitNormal[1];
            Jacobian_j[3][0] = Area*UnitNormal[2];
            Jacobian_j[4][0] = 0.0;

            Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitNormal[0] + (-(ru*rv)/(r*r))*Area*UnitNormal[1] + (-(ru*rw)/(r*r))*Area*UnitNormal[2];
            Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1] + (rw/r)*Area*UnitNormal[2];
            Jacobian_j[2][1] = (dp_drv)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[1];
            Jacobian_j[3][1] = (dp_drw)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[2];
            Jacobian_j[4][1] = (dp_drE)*Area*UnitNormal[0];

            Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitNormal[0] + (-(rv*rv)/(r*r) + dp_dr)*Area*UnitNormal[1] + (-(rv*rw)/(r*r))*Area*UnitNormal[2];
            Jacobian_j[1][2] = (rv/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[1];
            Jacobian_j[2][2] = (ru/r)*Area*UnitNormal[0] + (2*rv/r + dp_drv)*Area*UnitNormal[1] + (rw/r)*Area*UnitNormal[2];
            Jacobian_j[3][2] = (dp_drw)*Area*UnitNormal[1] + (rv/r)*Area*UnitNormal[2];
            Jacobian_j[4][2] = (dp_drE)*Area*UnitNormal[1];

            Jacobian_j[0][3] = (-(ru*rw)/(r*r))*Area*UnitNormal[0] + (-(rv*rw)/(r*r))*Area*UnitNormal[1] + (-(rw*rw)/(r*r) + dp_dr)*Area*UnitNormal[2];
            Jacobian_j[1][3] = (rw/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[2];
            Jacobian_j[2][3] = (rw/r)*Area*UnitNormal[1] + (dp_drv)*Area*UnitNormal[2];
            Jacobian_j[3][3] = (ru/r)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1] + (2*rw/r + dp_drw)*Area*UnitNormal[2];
            Jacobian_j[4][3] = (dp_drE)*Area*UnitNormal[2];

            Jacobian_j[0][4] = (ru*dH_dr)*Area*UnitNormal[0] + (rv*dH_dr)*Area*UnitNormal[1] + (rw*dH_dr)*Area*UnitNormal[2];
            Jacobian_j[1][4] = (H + ru*dH_dru)*Area*UnitNormal[0] + (rv*dH_dru)*Area*UnitNormal[1] + (rw*dH_dru)*Area*UnitNormal[2];
            Jacobian_j[2][4] = (ru*dH_drv)*Area*UnitNormal[0] + (H + rv*dH_drv)*Area*UnitNormal[1] + (rw*dH_drv)*Area*UnitNormal[2];
            Jacobian_j[3][4] = (ru*dH_drw)*Area*UnitNormal[0] + (rv*dH_drw)*Area*UnitNormal[1] + (H + rw*dH_drw)*Area*UnitNormal[2];
            Jacobian_j[4][4] = (ru*dH_drE)*Area*UnitNormal[0] + (rv*dH_drE)*Area*UnitNormal[1] + (rw*dH_drE)*Area*UnitNormal[2];
          }

          /*--- Mach number sensitivity ---*/

          USens[0] = 0.0; USens[1] = ru/Mach_Inf; USens[2] = rv/Mach_Inf;
          if (nDim == 2) { USens[3] = Gamma*Mach_Inf*p; }
          else { USens[3] = rw/Mach_Inf; USens[4] = Gamma*Mach_Inf*p; }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_Mach[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }

          /*--- AoA sensitivity ---*/

          USens[0] = 0.0;
          if (nDim == 2) { USens[1] = -rv; USens[2] = ru; USens[3] = 0.0; }
          else { USens[1] = -rw; USens[2] = 0.0; USens[3] = ru; USens[4] = 0.0; }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_AoA[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }

          /*--- Pressure sensitivity ---*/

          USens[0] = r/p; USens[1] = ru/p; USens[2] = rv/p;
          if (nDim == 2) { USens[3] = rE/p; }
          else { USens[3] = rw/p; USens[4] = rE/p; }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_Press[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }

          /*--- Temperature sensitivity ---*/

          T = p/(r*Gas_Constant);
          USens[0] = -r/T; USens[1] = 0.5*ru/T; USens[2] = 0.5*rv/T;
          if (nDim == 2) { USens[3] = (ru*ru + rv*rv + rw*rw)/(r*T); }
          else { USens[3] = 0.5*rw/T; USens[4] = (ru*ru + rv*rv + rw*rw)/(r*T); }
          for (iPos = 0; iPos < nVar; iPos++) {
            for (jPos = 0; jPos < nVar; jPos++) {
              Sens_Temp[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
            }
          }
        }
      }

      Total_Sens_Mach  -= Sens_Mach[iMarker] * scale * factor;
      Total_Sens_AoA   -= Sens_AoA[iMarker] * scale * factor;
      Total_Sens_Press -= Sens_Press[iMarker] * scale * factor;
      Total_Sens_Temp  -= Sens_Temp[iMarker] * scale * factor;

    }

  }

  /*--- Explicit contribution from objective function quantity ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)) {

      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      Sens_BPress[iMarker] = 0.0;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          p = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          d = nodes->GetForceProj_Vector(iPoint);
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

          /*--- Mach number sensitivity ---*/

          for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = -(2.0/Mach_Inf)*d[iPos];
          for (iPos = 0; iPos < nDim; iPos++) Sens_Mach[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

          /*--- AoA sensitivity ---*/

          if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT ||
              config->GetKind_ObjFunc() == LIFT_COEFFICIENT ||
              config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT ||
              config->GetKind_ObjFunc() == EQUIVALENT_AREA ||
              config->GetKind_ObjFunc() == NEARFIELD_PRESSURE) {
            if (nDim == 2) {
              D[0][0] = 0.0; D[0][1] = -1.0;
              D[1][0] = 1.0; D[1][1] = 0.0;
            }
            else {
              D[0][0] = 0.0; D[0][1] = 0.0; D[0][2] = -1.0;
              D[1][0] = 0.0; D[1][1] = 0.0; D[1][2] = 0.0;
              D[2][0] = 1.0; D[2][1] = 0.0; D[2][2] = 0.0;
            }
            for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = 0.0;
            for (iPos = 0; iPos < nDim; iPos++) {
              for (jPos = 0; jPos < nDim; jPos++)
                Dd[iPos] += D[iPos][jPos]*d[jPos];
            }
          }

          /*--- Coefficients with no explicit AoA dependece ---*/

          else {
            for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
          }

          for (iPos = 0; iPos < nDim; iPos++)
            Sens_AoA[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

          /*--- Pressure sensitivity ---*/

          for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = -(1/p)*d[iPos];
          for (iPos = 0; iPos<nDim; iPos++)
            Sens_Press[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

          /*--- Temperature sensitivity ---*/

          for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
          for (iPos = 0; iPos<nDim; iPos++)
            Sens_Temp[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];

        }
      }

      Total_Sens_Mach  += Sens_Mach[iMarker] * scale * factor;
      Total_Sens_AoA   += Sens_AoA[iMarker] * scale * factor;
      Total_Sens_Press += Sens_Press[iMarker] * scale * factor;
      Total_Sens_Temp  += Sens_Temp[iMarker] * scale * factor;

    }
  }

  
#ifdef HAVE_MPI
  
  su2double MyTotal_Sens_Geo   = Total_Sens_Geo;     Total_Sens_Geo = 0.0;
  su2double MyTotal_Sens_Mach  = Total_Sens_Mach;    Total_Sens_Mach = 0.0;
  su2double MyTotal_Sens_AoA   = Total_Sens_AoA;     Total_Sens_AoA = 0.0;
  su2double MyTotal_Sens_Press = Total_Sens_Press;   Total_Sens_Press = 0.0;
  su2double MyTotal_Sens_Temp  = Total_Sens_Temp;    Total_Sens_Temp = 0.0;
  
  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif

  delete [] USens;
  delete [] UnitNormal;
  delete [] normal_grad_vel;
  delete [] tang_deriv_psi5;
  delete [] tang_deriv_T;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Sigma[iDim];
  delete [] Sigma;
  delete [] normal_grad_gridvel;
  delete [] normal_grad_v_ux;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Sigma_Psi5v[iDim];
  delete [] Sigma_Psi5v;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] tau[iDim];
  delete [] tau;
  delete [] Velocity;

}

void CAdjNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar, jDim;
  unsigned long iVertex, iPoint, total_index, Point_Normal;
  
  su2double *d, l1psi, vartheta, Sigma_5, **PsiVar_Grad, phi[3] = {};
  su2double sq_vel, ProjGridVel, Enthalpy = 0.0, *GridVel;
  su2double ViscDens, XiDens, Density, Pressure = 0.0, dPhiE_dn;
  su2double Laminar_Viscosity = 0.0, Eddy_Viscosity = 0.0;
  su2double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz;
  su2double Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5;
  su2double Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  su2double *Coord_i, *Coord_j, dist_ij_2;
  su2double dSigmaxx_phi1, dSigmayy_phi1, dSigmazz_phi1, dSigmaxy_phi1, dSigmaxz_phi1, dSigmayz_phi1;
  su2double dSigmaxx_phi2, dSigmayy_phi2, dSigmazz_phi2, dSigmaxy_phi2, dSigmaxz_phi2, dSigmayz_phi2;
  su2double dSigmaxx_phi3, dSigmayy_phi3, dSigmazz_phi3, dSigmaxy_phi3, dSigmaxz_phi3, dSigmayz_phi3;
  su2double dSigma5_psi5;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  
  su2double *Psi = new su2double[nVar];
  su2double **Tau = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Tau[iDim] = new su2double [nDim];
  su2double *Velocity = new su2double[nDim];
  su2double *Normal = new su2double[nDim];
  su2double *Edge_Vector = new su2double[nDim];
  su2double **GradPhi = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new su2double [nDim];
  su2double *GradPsiE = new su2double [nDim];
  
  /*--- Loop over all of the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0; Res_Visc_i[iVar] = 0.0;
        if (implicit) { for (jVar = 0; jVar < nVar; jVar ++) Jacobian_ii[iVar][jVar] = 0.0; }
      }
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Psi[iVar] = nodes->GetSolution(iPoint,iVar);
      
      /*--- Get the force projection vector (based on the objective function) ---*/
      
      d = nodes->GetForceProj_Vector(iPoint);
      
      /*--- Set the adjoint velocity BC ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) { phi[iDim] = d[iDim]; }
      
      /*--- Correct the adjoint velocity BC for dynamic meshes ---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          phi[iDim] -= Psi[nDim+1]*GridVel[iDim];
      }
      
      /*--- Impose the value of the adjoint velocity as a strong boundary
       condition (Dirichlet). Fix the adjoint velocity and remove any addtional
       contribution to the residual at this node. ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        nodes->SetSolution_Old(iPoint,iDim+1, phi[iDim]);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);
      
      /*--- Compute additional contributions to the adjoint density and energy
       equations which will be added to the residual (weak imposition) ---*/

      /*--- Energy residual due to the convective term ---*/

      l1psi = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        l1psi += Normal[iDim]*d[iDim];
      Res_Conv_i[nDim+1] = l1psi*Gamma_Minus_One;

      /*--- Flux contribution and Jacobian contributions for moving
         walls. Note that these are only for the adjoint density and
         adjoint energy equations (the adjoint vel. uses a strong BC). ---*/

      if (grid_movement) {

        /*--- Get the grid velocity at this node and impose v = u_wall ---*/

        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Velocity[iDim] = GridVel[iDim];

        /*--- Get some additional quantities from the flow solution ---*/

        Density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        Pressure = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);
        Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);
        Laminar_Viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        Eddy_Viscosity = solver_container[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint); // Should be zero at the wall

        ViscDens = (Laminar_Viscosity + Eddy_Viscosity) / Density;
        XiDens = Gamma * (Laminar_Viscosity/Prandtl_Lam + Eddy_Viscosity/Prandtl_Turb) / Density;

        /*--- Compute projections, velocity squared divided by two, and
           other inner products. Note that we are imposing v = u_wall from
           the direct problem and that phi = d - \psi_5 * v ---*/

        ProjGridVel = 0.0; sq_vel = 0.0;
        vartheta = Psi[0] + Psi[nDim+1]*Enthalpy;
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjGridVel += GridVel[iDim]*Normal[iDim];
          sq_vel      += 0.5*GridVel[iDim]*GridVel[iDim];
          vartheta    += GridVel[iDim]*phi[iDim];
        }

        /*--- Convective flux at the wall node (adjoint density) ---*/

        Res_Conv_i[0] = -vartheta*ProjGridVel + l1psi*Gamma_Minus_One*sq_vel;

        /*--- Implicit contributions from convective part ---*/

        if (implicit) {
          Jacobian_ii[0][0] += -ProjGridVel;
          Jacobian_ii[0][nVar-1] += -ProjGridVel * Enthalpy;
        }

        /*--- Viscous flux contributions at the wall node. Impose dPhiE_dn = 0
           (adiabatic walls with frozen viscosity). ---*/

        dPhiE_dn = 0.0;

        /*--- Store the adjoint velocity and energy gradients for clarity ---*/

        PsiVar_Grad = nodes->GetGradient(iPoint);
        for (iDim = 0; iDim < nDim; iDim++) {
          GradPsiE[iDim] =  PsiVar_Grad[nVar-1][iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            GradPhi[iDim][jDim] =  PsiVar_Grad[iDim+1][jDim];
        }

        if (nDim == 2) {

          /*--- Compute the adjoint stress tensor ---*/

          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5;
          eta_yy    = Sigma_yy + Sigma_yy5;
          eta_xy    = Sigma_xy + Sigma_xy5;

          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/

          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy
              + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
              - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[nDim+1] = Sigma_5;

          /*--- Computation of the Jacobians at Point i---*/

          if (implicit) {

            /*--- Compute closest normal neighbor ---*/

            Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

            /*--- Get coordinates of i & nearest normal and compute distance ---*/

            Coord_i = geometry->node[iPoint]->GetCoord();
            Coord_j = geometry->node[Point_Normal]->GetCoord();
            dist_ij_2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
              dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
            }

            dSigmaxx_phi1 = -FOUR3 * ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmaxx_phi2 =   TWO3 * ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmayy_phi1 =   TWO3 * ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmayy_phi2 = -FOUR3 * ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmaxy_phi1 = -ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmaxy_phi2 = -ViscDens * Edge_Vector[0]/dist_ij_2;

            //              dSigmaxx5_psi5 = -ViscDens * ( FOUR3*Velocity[0]*Edge_Vector[0] -  TWO3*Velocity[1]*Edge_Vector[1] )/dist_ij_2;
            //              dSigmayy5_psi5 = -ViscDens * (- TWO3*Velocity[0]*Edge_Vector[0] + FOUR3*Velocity[1]*Edge_Vector[1] )/dist_ij_2;
            //              dSigmaxy5_psi5 = -ViscDens * ( Velocity[0]*Edge_Vector[1] + Velocity[1]*Edge_Vector[0] )/dist_ij_2;
            dSigma5_psi5   = -XiDens * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;

            Jacobian_ii[0][0] += 0.0;
            Jacobian_ii[0][1] += -( Velocity[0]*Normal[0]*dSigmaxx_phi1 + Velocity[1]*Normal[1]*dSigmayy_phi1
                + (Velocity[0]*Normal[1] + Velocity[1]*Normal[0])*dSigmaxy_phi1 );
            Jacobian_ii[0][2] += -( Velocity[0]*Normal[0]*dSigmaxx_phi2 + Velocity[1]*Normal[1]*dSigmayy_phi2
                + (Velocity[0]*Normal[1] + Velocity[1]*Normal[0])*dSigmaxy_phi2 );
            Jacobian_ii[0][3] += (sq_vel - Pressure/(Density*Gamma_Minus_One)) * dSigma5_psi5;

            Jacobian_ii[3][0] += 0.0;
            Jacobian_ii[3][1] += 0.0;
            Jacobian_ii[3][2] += 0.0;
            Jacobian_ii[3][3] += dSigma5_psi5;

          }


        } else if (nDim == 3) {

          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_zz  = ViscDens * (-TWO3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] + FOUR3 * GradPhi[2][2]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xz  = ViscDens * (GradPhi[2][0] + GradPhi[0][2]);
          Sigma_yz  = ViscDens * (GradPhi[2][1] + GradPhi[1][2]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_zz5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] + FOUR3 * Velocity[2] * GradPsiE[2]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_xz5 = ViscDens * (Velocity[0] * GradPsiE[2] + Velocity[2] * GradPsiE[0]);
          Sigma_yz5 = ViscDens * (Velocity[1] * GradPsiE[2] + Velocity[2] * GradPsiE[1]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
          eta_xy    = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/

          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy + Velocity[2] * Normal[2] * eta_zz
              + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
              + (Velocity[0] * Normal[2] + Velocity[2] * Normal[0]) * eta_xz
              + (Velocity[2] * Normal[1] + Velocity[1] * Normal[2]) * eta_yz
              - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[nDim+1] = Sigma_5;

          /*--- Computation of the Jacobians at Point i---*/

          if (implicit) {

            /*--- Compute closest normal neighbor ---*/

            Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

            /*--- Get coordinates of i & nearest normal and compute distance ---*/

            Coord_i = geometry->node[iPoint]->GetCoord();
            Coord_j = geometry->node[Point_Normal]->GetCoord();
            dist_ij_2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
              dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
            }

            dSigmaxx_phi1 = -FOUR3 * ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmaxx_phi2 =   TWO3 * ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmaxx_phi3 =   TWO3 * ViscDens * Edge_Vector[2]/dist_ij_2;
            dSigmayy_phi1 =   TWO3 * ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmayy_phi2 = -FOUR3 * ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmayy_phi3 =   TWO3 * ViscDens * Edge_Vector[2]/dist_ij_2;
            dSigmazz_phi1 =   TWO3 * ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmazz_phi2 =   TWO3 * ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmazz_phi3 = -FOUR3 * ViscDens * Edge_Vector[2]/dist_ij_2;
            dSigmaxy_phi1 = -ViscDens * Edge_Vector[1]/dist_ij_2;
            dSigmaxy_phi2 = -ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmaxy_phi3 = 0;
            dSigmaxz_phi1 = -ViscDens * Edge_Vector[2]/dist_ij_2;
            dSigmaxz_phi2 = 0;
            dSigmaxz_phi3 = -ViscDens * Edge_Vector[0]/dist_ij_2;
            dSigmayz_phi1 = 0;
            dSigmayz_phi2 = -ViscDens * Edge_Vector[2]/dist_ij_2;
            dSigmayz_phi3 = -ViscDens * Edge_Vector[1]/dist_ij_2;

            //              dSigmaxx5_psi5 = -ViscDens * ( FOUR3*Velocity[0]*Edge_Vector[0] -  TWO3*Velocity[1]*Edge_Vector[1] -  TWO3*Velocity[2]*Edge_Vector[2])/dist_ij_2;
            //              dSigmayy5_psi5 = -ViscDens * (- TWO3*Velocity[0]*Edge_Vector[0] + FOUR3*Velocity[1]*Edge_Vector[1] -  TWO3*Velocity[2]*Edge_Vector[2])/dist_ij_2;
            //              dSigmazz5_psi5 = -ViscDens * (- TWO3*Velocity[0]*Edge_Vector[0] -  TWO3*Velocity[1]*Edge_Vector[1] + FOUR3*Velocity[2]*Edge_Vector[2])/dist_ij_2;
            //              dSigmaxy5_psi5 = -ViscDens * ( Velocity[0]*Edge_Vector[1] + Velocity[1]*Edge_Vector[0] )/dist_ij_2;
            //              dSigmaxz5_psi5 = -ViscDens * ( Velocity[0]*Edge_Vector[2] + Velocity[2]*Edge_Vector[0] )/dist_ij_2;
            //              dSigmayz5_psi5 = -ViscDens * ( Velocity[1]*Edge_Vector[2] + Velocity[2]*Edge_Vector[1] )/dist_ij_2;
            dSigma5_psi5   = -XiDens * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;

            Jacobian_ii[0][0] += 0.0;
            Jacobian_ii[0][1] += -( Velocity[0]*Normal[0]*dSigmaxx_phi1 + Velocity[1]*Normal[1]*dSigmayy_phi1 + Velocity[2]*Normal[2]*dSigmazz_phi1
                + (Velocity[0]*Normal[1] + Velocity[1]*Normal[0])*dSigmaxy_phi1
                + (Velocity[0]*Normal[2] + Velocity[2]*Normal[0])*dSigmaxz_phi1
                + (Velocity[2]*Normal[1] + Velocity[1]*Normal[2])*dSigmayz_phi1 );
            Jacobian_ii[0][2] += -( Velocity[0]*Normal[0]*dSigmaxx_phi2 + Velocity[1]*Normal[1]*dSigmayy_phi2 + Velocity[2]*Normal[2]*dSigmazz_phi2
                + (Velocity[0]*Normal[1] + Velocity[1]*Normal[0])*dSigmaxy_phi2
                + (Velocity[0]*Normal[2] + Velocity[2]*Normal[0])*dSigmaxz_phi2
                + (Velocity[2]*Normal[1] + Velocity[1]*Normal[2])*dSigmayz_phi2 );
            Jacobian_ii[0][3] += -( Velocity[0]*Normal[0]*dSigmaxx_phi3 + Velocity[1]*Normal[1]*dSigmayy_phi3 + Velocity[2]*Normal[2]*dSigmazz_phi3
                + (Velocity[0]*Normal[1] + Velocity[1]*Normal[0])*dSigmaxy_phi3
                + (Velocity[0]*Normal[2] + Velocity[2]*Normal[0])*dSigmaxz_phi3
                + (Velocity[2]*Normal[1] + Velocity[1]*Normal[2])*dSigmayz_phi3 );
            Jacobian_ii[0][4] += (sq_vel - Pressure/(Density*Gamma_Minus_One)) * dSigma5_psi5;

            Jacobian_ii[4][0] += 0.0;
            Jacobian_ii[4][1] += 0.0;
            Jacobian_ii[4][2] += 0.0;
            Jacobian_ii[4][3] += 0.0;
            Jacobian_ii[4][4] += dSigma5_psi5;

          }
        }
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
    
  }
  
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Tau[iDim];
  delete [] Tau;
  delete [] Psi;
  delete [] Velocity;
  delete [] Normal;
  delete [] Edge_Vector;
  delete [] GradPsiE;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] GradPhi[iDim];
  delete [] GradPhi;
  
}


void CAdjNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, total_index;
  unsigned short iDim, iVar, jVar, jDim;
  su2double *d, q, *U, dVisc_T, rho, pressure, div_phi,
  force_stress, Sigma_5, **PsiVar_Grad, phi[3] = {0.0,0.0,0.0};
  su2double phis1, phis2, sq_vel, ProjVel, Enthalpy, *GridVel, phi_u, d_n;
  su2double Energy, ViscDens, XiDens, Density, SoundSpeed, Pressure, dPhiE_dn, Laminar_Viscosity, Eddy_Viscosity,
  Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  su2double kGTdotn=0.0, Area=0.0, Xi=0.0;
  
  su2double *Psi = new su2double[nVar];
  su2double **Tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Tau[iDim] = new su2double [nDim];
  su2double *Velocity = new su2double[nDim];
  su2double *Normal = new su2double[nDim];
  
  su2double **GradPhi = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new su2double [nDim];
  su2double *GradPsiE = new su2double [nDim];
  su2double *GradT;// = new su2double[nDim];
  su2double *GradP;
  su2double *GradDens;
  su2double *dPoRho2 = new su2double[nDim];
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool heat_flux_obj;
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Thermal_Conductivity;
  su2double invrho3;
  su2double Volume;
  su2double mu2;
  su2double gpsiAv2;
  su2double gpsi5n;

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  string Monitoring_Tag;
  unsigned short jMarker, iMarker_Monitoring=0;
  su2double Weight_ObjFunc = 1.0;
  
  /*--- Identify marker monitoring index ---*/
  for (jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
    Monitoring_Tag = config->GetMarker_Monitoring_TagBound(jMarker);
    if (Monitoring_Tag==Marker_Tag)
      iMarker_Monitoring = jMarker;
  }
  /*-- Get objective weight --*/
  Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
  heat_flux_obj  = ((config->GetKind_ObjFunc(iMarker_Monitoring) == TOTAL_HEATFLUX) ||
                           (config->GetKind_ObjFunc(iMarker_Monitoring) == MAXIMUM_HEATFLUX) ||
                           (config->GetKind_ObjFunc(iMarker_Monitoring) == INVERSE_DESIGN_HEATFLUX));

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0;
        Res_Visc_i[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_ii[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Psi[iVar] = nodes->GetSolution(iPoint,iVar);
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      Volume = geometry->node[iPoint]->GetVolume();
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Get the force projection vector (based on the objective function) ---*/
      d = nodes->GetForceProj_Vector(iPoint);
      
      /*--- Adjustments to strong boundary condition for dynamic meshes ---*/
      if ( grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim] - Psi[nVar-1]*GridVel[iDim];
        }
      } else {
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim];
        }
      }
      
      /*--- Strong BC imposition for the adjoint velocity equations ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);
      for (iDim = 0; iDim < nDim; iDim++)
        nodes->SetSolution_Old(iPoint,iDim+1, phi[iDim]);
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
      /*--- Get transport coefficient information ---*/
      Laminar_Viscosity    = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
      Eddy_Viscosity       = solver_container[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint);
      Thermal_Conductivity = Cp * ( Laminar_Viscosity/Prandtl_Lam
                                   +Eddy_Viscosity/Prandtl_Turb);
      
//      GradV = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
      
      /*--- Calculate Dirichlet condition for energy equation ---*/
      if (!heat_flux_obj) {
        q = 0.0;
      }
      else {

        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);

        /*--- Temperature gradient term ---*/
        GradT = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint)[0];
        kGTdotn = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          kGTdotn += Cp * Laminar_Viscosity/Prandtl_Lam*GradT[iDim]*Normal[iDim]/Area;
        // Cp * Viscosity/Prandtl_Lam matches term used in solver_direct_mean
        /*--- constant term to multiply max heat flux objective ---*/
        Xi = solver_container[FLOW_SOL]->GetTotal_HeatFlux(); // versions for max heat flux
        Xi = pow(Xi, 1.0/pnorm-1.0)/pnorm;

        /*--- Boundary condition value ---*/
        q = Xi * pnorm * pow(kGTdotn, pnorm-1.0)*Area*Weight_ObjFunc;
      }
      
      /*--- Strong BC enforcement of the energy equation ---*/
      LinSysRes.SetBlock_Zero(iPoint, nVar-1);
      nodes->SetEnergy_ResTruncError_Zero(iPoint);
      nodes->SetSolution_Old(iPoint,nDim+1, q);
      if (implicit) {
        iVar = nDim+1;
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
      /*--- Additional contributions to adjoint density (weak imposition) ---*/

      /*--- Acquire gradient information ---*/
      PsiVar_Grad = nodes->GetGradient(iPoint);
      GradP    = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint)[nVar-1];
      GradDens = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint)[nVar];

      /*--- Acqure flow information ---*/
      rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
      pressure = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);
      invrho3 = (1.0/rho)*(1.0/rho)*(1.0/rho);

      /*--- Calculate supporting quantities ---*/
      mu2 = Thermal_Conductivity/Cp;
      gpsiAv2 = 0.0;
      gpsi5n = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dPoRho2[iDim] = (GradP[iDim]*rho - 2.0*GradDens[iDim]*pressure)*invrho3;
        gpsiAv2 += -mu2*Gamma/Gamma_Minus_One * PsiVar_Grad[nVar-1][iDim]*dPoRho2[iDim];
        gpsi5n += PsiVar_Grad[nVar-1][iDim]*Normal[iDim];
      }

      /*--- Apply first order term to boundary ---*/
      Res_Conv_i[0] = gpsiAv2*Volume;

      /*--- Apply second order term to boundary ---*/
      Res_Visc_i[0] = -mu2*Gamma/(rho*Gamma_Minus_One)*(pressure/rho)*gpsi5n;

      /*--- Components of the effective and adjoint stress tensors ---*/
      PsiVar_Grad = nodes->GetGradient(iPoint);
      div_phi = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        div_phi += PsiVar_Grad[iDim+1][iDim];
        for (jDim = 0; jDim < nDim; jDim++)
          Tau[iDim][jDim] = (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
      }
      for (iDim = 0; iDim < nDim; iDim++)
        Tau[iDim][iDim] -= TWO3*div_phi;

      /*--- force_stress = n_i \Tau_{ij} d_j ---*/
      force_stress = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          force_stress += Normal[iDim]*Tau[iDim][jDim]*d[jDim];

      /*--- \partial \mu_dyn \partial T ---*/
      //        mu_dyn = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
      //        Temp = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
      dVisc_T = 0.0;  // dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));

      /*--- \Sigma_5 Check Area computation for Res_Conv[0] ---*/
      Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;

      /*--- Imposition of residuals ---*/
      rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
      pressure = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);
      Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);

      /*--- Flux contribution and Jacobian contributions for moving
         walls. Note that these are only for the adjoint density and
         adjoint energy equations (the adjoint vel. uses a strong BC). ---*/
      if (grid_movement) {

        /*--- Get the appropriate grid velocity at this node ---*/
        GridVel = geometry->node[iPoint]->GetGridVel();

        /*--- Get the enthalpy from the direct solution ---*/
        Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);

        /*--- Compute projections, velocity squared divided by two, and
           other inner products. Note that we are imposing v = u_wall from
           the direct problem and that phi = d - \psi_5 * v ---*/
        ProjVel = 0.0; sq_vel = 0.0; phi_u = 0.0; d_n = 0.0;
        phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjVel += GridVel[iDim]*Normal[iDim];
          sq_vel  += 0.5*GridVel[iDim]*GridVel[iDim];
          phis1   += Normal[iDim]*phi[iDim];
          phis2   += GridVel[iDim]*phi[iDim];
          phi_u   += GridVel[iDim]*phi[iDim];
          d_n     += d[iDim]*Normal[iDim];
        }
        //          phis1 += ProjVel * Psi[nVar-1];

        /*--- Convective flux at the wall node (adjoint density & energy only) ---*/

        /*--- Version 1 (full) ---*/
        //Res_Conv_i[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel - ProjVel*Psi[0];
        //Res_Conv_i[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One - ProjVel*Psi[nVar-1];

        /*--- Simplified version ---*/
        Res_Conv_i[0] = -(Psi[0] + phi_u + Psi[nVar-1]*Enthalpy)*ProjVel + d_n*Gamma_Minus_One*sq_vel;

        /*--- TO DO: Implicit contributions for convective part ---*/


        /*--- Viscous flux contributions at the wall node ---*/
        U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
        Laminar_Viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        Eddy_Viscosity = solver_container[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint); // Should be zero at the wall
        Density = U[0];
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = GridVel[iDim];
        }
        Energy = U[nDim+1] / Density;
        SoundSpeed = sqrt(Gamma*Gamma_Minus_One*(Energy-sq_vel));
        Pressure = (SoundSpeed * SoundSpeed * Density) / Gamma;
        ViscDens = (Laminar_Viscosity + Eddy_Viscosity) / Density;
        XiDens = Gamma * (Laminar_Viscosity/Prandtl_Lam + Eddy_Viscosity/Prandtl_Turb) / Density;

        /*--- Average of the derivatives of the adjoint variables ---*/
        PsiVar_Grad = nodes->GetGradient(iPoint);

        for (iDim = 0; iDim < nDim; iDim++) {
          GradPsiE[iDim] =  PsiVar_Grad[nVar-1][iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            GradPhi[iDim][jDim] =  PsiVar_Grad[iDim+1][jDim];
        }

        /*--- Impose dPhiE_dn = 0 (adiabatic walls with frozen viscosity). Note
           that this is where a different adjoint boundary condition for temperature
           could be imposed. ---*/
        dPhiE_dn = 0.0;

        if (nDim ==2) {

          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5;
          eta_yy    = Sigma_yy + Sigma_yy5;
          eta_xy    = Sigma_xy + Sigma_xy5;

          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy
              + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
              - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[1] = 0.0;
          Res_Visc_i[2] = 0.0;

        } else if (nDim == 3) {

          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_zz  = ViscDens * (-TWO3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] + FOUR3 * GradPhi[2][2]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xz  = ViscDens * (GradPhi[2][0] + GradPhi[0][2]);
          Sigma_yz  = ViscDens * (GradPhi[2][1] + GradPhi[1][2]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_zz5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] + FOUR3 * Velocity[2] * GradPsiE[2]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_xz5 = ViscDens * (Velocity[0] * GradPsiE[2] + Velocity[2] * GradPsiE[0]);
          Sigma_yz5 = ViscDens * (Velocity[1] * GradPsiE[2] + Velocity[2] * GradPsiE[1]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
          eta_xy    = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy + Velocity[2] * Normal[2] * eta_zz
              + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
              + (Velocity[0] * Normal[2] + Velocity[2] * Normal[0]) * eta_xz
              + (Velocity[2] * Normal[1] + Velocity[1] * Normal[2]) * eta_yz
              - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[1] = 0.0;
          Res_Visc_i[2] = 0.0;
          Res_Visc_i[3] = 0.0;
        }
      }
      
      /*--- Update convective and viscous residuals ---*/
      LinSysRes.AddBlock(iPoint, Res_Conv_i);
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      }
      
    }
    
  }
  
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Tau[iDim];
  delete [] Tau;
  delete [] Psi;
  delete [] Velocity;
  delete [] Normal;
  delete [] GradPsiE;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] GradPhi[iDim];
  delete [] GradPhi;
  delete [] dPoRho2;
}

