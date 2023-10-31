/*!
 * \file CAdjEulerSolver.cpp
 * \brief Main subroutines for solving Euler adjoint problems.
 * \author F. Palacios, T. Economon, H. Kline
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


#include "../../include/solvers/CAdjEulerSolver.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CAdjEulerSolver::CAdjEulerSolver() : CSolver() {

  /*--- Array initialization ---*/
  Phi_Inf = nullptr;
  Sens_Mach = nullptr;
  Sens_AoA = nullptr;
  Sens_Geo = nullptr;
  Sens_Press = nullptr;
  Sens_Temp = nullptr;
  Sens_BPress = nullptr;
  Jacobian_Axisymmetric = nullptr;
  CSensitivity = nullptr;
  FlowPrimVar_i = nullptr;
  FlowPrimVar_j = nullptr;
  DonorAdjVar = nullptr;
  DonorGlobalIndex = nullptr;

}

CAdjEulerSolver::CAdjEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  unsigned long iPoint, iVertex, iMarker, jMarker;
  string text_line, mesh_filename;
  unsigned short iDim, iVar;
  ifstream restart_file;
  string filename, AdjExt;
  su2double myArea_Monitored, *Normal;

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
  unsigned short iMarker_Monitoring, ObjFunc;
  bool grid_movement  = config->GetGrid_Movement();

  /*--- Array initialization ---*/

  Phi_Inf = nullptr;
  Sens_Mach = nullptr;
  Sens_AoA = nullptr;
  Sens_Geo = nullptr;
  Sens_Press = nullptr;
  Sens_Temp = nullptr;
  Sens_BPress = nullptr;
  Jacobian_Axisymmetric = nullptr;
  CSensitivity = nullptr;
  FlowPrimVar_i = nullptr;
  FlowPrimVar_j = nullptr;
  DonorAdjVar = nullptr;
  DonorGlobalIndex = nullptr;

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

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
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
  FlowPrimVar_i = new su2double[nDim+7]();
  FlowPrimVar_j = new su2double[nDim+7]();

  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) {
    iPoint_UndLapl.resize(nPoint);
    jPoint_UndLapl.resize(nPoint);
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

    if (axisymmetric) {
      Jacobian_Axisymmetric = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_Axisymmetric[iVar] = new su2double [nVar];
    }
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
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

  space_centered = config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED;


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
            if (geometry->nodes->GetDomain(iPoint)) {
              Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
              myArea_Monitored += GeometryToolbox::Norm(nDim, Normal);
            }
          }
          break;
        }
      }
    }
  }

#ifdef HAVE_MPI
  Area_Monitored = 0.0;
  SU2_MPI::Allreduce(&myArea_Monitored, &Area_Monitored, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
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

  SolverName = "ADJ.FLOW";
}

CAdjEulerSolver::~CAdjEulerSolver() {
  unsigned short iVar, iMarker;

  delete [] Phi_Inf;
  delete [] Sens_Mach;
  delete [] Sens_AoA;
  delete [] Sens_Geo;
  delete [] Sens_Press;
  delete [] Sens_Temp;
  delete [] Sens_BPress;
  delete [] FlowPrimVar_i;
  delete [] FlowPrimVar_j;

  if (Jacobian_Axisymmetric != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_Axisymmetric[iVar];
    delete [] Jacobian_Axisymmetric;
  }

  if (CSensitivity != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CSensitivity[iMarker];
    delete [] CSensitivity;
  }

  delete nodes;
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

  /*--- MPI status and request arrays for non-blocking communications ---*/

  SU2_MPI::Status status;
  SU2_MPI::Request req;

  /*--- Define buffer vector interior domain ---*/

  su2double        *Buffer_Send_AdjVar          = nullptr;
  auto        *iAdjVar          = new su2double [nVar];

  auto *nPointTotal_s = new unsigned long[size];
  auto *nPointTotal_r = new unsigned long[size];

  unsigned long Buffer_Size_AdjVar          = 0;
  unsigned long PointTotal_Counter = 0;

  /*--- Allocate the memory that we only need if we have MPI support ---*/

  su2double        *Buffer_Receive_AdjVar          = nullptr;

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
          if ((iDomain == jDomain) && (geometry->nodes->GetDomain(iPoint))) {
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

      /*--- Communicate the counts to iDomain with non-blocking sends ---*/

      SU2_MPI::Isend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, SU2_MPI::GetComm(), &req);
      SU2_MPI::Request_free(&req);

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

          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/

          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, SU2_MPI::GetComm(), &status);

        }
      }

    }
  }

  /*--- Wait for the non-blocking sends to complete. ---*/

  SU2_MPI::Barrier(SU2_MPI::GetComm());

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
          if ((iDomain == jDomain) && (geometry->nodes->GetDomain(iPoint))) {

            iGlobalIndex = geometry->nodes->GetGlobalIndex(iPoint);
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

      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/

      SU2_MPI::Isend(&Buffer_Send_AdjVar[PointTotal_Counter*(nVar+3)],
                     nPointTotal_s[iDomain]*(nVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  SU2_MPI::GetComm(), &req);
      SU2_MPI::Request_free(&req);
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

  SU2_MPI::Barrier(SU2_MPI::GetComm());

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
                    iDomain, rank, SU2_MPI::GetComm(), &status);

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

  SU2_MPI::Barrier(SU2_MPI::GetComm());

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
  const auto RefOriginMoment = config->GetRefOriginMoment(0);
  su2double dCD_dCL          = config->GetdCD_dCL();
  su2double dCMx_dCL         = config->GetdCMx_dCL();
  su2double dCMy_dCL         = config->GetdCMy_dCL();
  su2double dCMz_dCL         = config->GetdCMz_dCL();
  bool Fixed_CL              = config->GetFixed_CL_Mode();

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

        x = geometry->nodes->GetCoord(iPoint, 0);
        y = geometry->nodes->GetCoord(iPoint, 1);
        if (nDim == 3) z = geometry->nodes->GetCoord(iPoint, 2);

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

            }
            if (nDim == 3) {

              ForceProj_Vector[0] += Weight_ObjFunc*cos(Alpha)*cos(Beta);
              ForceProj_Vector[1] += Weight_ObjFunc*sin(Beta);
              ForceProj_Vector[2] += Weight_ObjFunc*sin(Alpha)*cos(Beta);

              /*--- Modification to run at a fixed CL value ---*/

              if (Fixed_CL) { ForceProj_Vector[0] += dCD_dCL*Weight_ObjFunc*sin(Alpha); ForceProj_Vector[1] -= 0.0; ForceProj_Vector[2] -= dCD_dCL*Weight_ObjFunc*cos(Alpha); }

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

void CAdjEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  const bool restart = config->GetRestart();
  const bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                          (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  if (restart) {
    for (auto iMesh = 1ul; iMesh <= config->GetnMGLevels(); iMesh++) {
      MultigridRestriction(*geometry[iMesh - 1], solver_container[iMesh - 1][ADJFLOW_SOL]->GetNodes()->GetSolution(),
                           *geometry[iMesh], solver_container[iMesh][ADJFLOW_SOL]->GetNodes()->GetSolution());
      solver_container[iMesh][ADJFLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][ADJFLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    }
  }

  /*--- The value of the solution for the first iteration of the dual time ---*/
  for (auto iMesh = 0ul; iMesh <= config->GetnMGLevels(); iMesh++) {
    if ((TimeIter == 0) && (dual_time)) {
      solver_container[iMesh][ADJFLOW_SOL]->GetNodes()->Set_Solution_time_n();
      solver_container[iMesh][ADJFLOW_SOL]->GetNodes()->Set_Solution_time_n1();
    }
  }

}

void CAdjEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  su2double SharpEdge_Distance;
  bool physical = true;

  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/

  bool implicit       = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool muscl          = config->GetMUSCL_AdjFlow();
  bool limiter        = (config->GetKind_SlopeLimit_AdjFlow() != LIMITER::NONE);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == CENTERED::JST);
  bool fixed_cl       = config->GetFixed_CL_Mode();
  bool eval_dof_dcx   = config->GetEval_dOF_dCX();

  /*--- Update the objective function coefficient to guarantee zero gradient. ---*/

  if (fixed_cl && eval_dof_dcx) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

  /*--- Residual initialization ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Get the distance form a sharp edge ---*/

    SharpEdge_Distance = geometry->nodes->GetSharpEdge_Distance(iPoint);

    /*--- Set the primitive variables compressible
     adjoint variables ---*/

    physical = nodes->SetPrimVar(iPoint,SharpEdge_Distance, false, config);

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;

    /*--- Initialize the convective residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  if ((muscl) && (iMesh == MESH_0)) {

    /*--- Gradient computation for MUSCL reconstruction. ---*/

    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetSolution_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);

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
    unsigned long MyErrorCounter = nonPhysicalPoints; nonPhysicalPoints = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &nonPhysicalPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(nonPhysicalPoints);
  }

}

void CAdjEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                        CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  unsigned long iEdge, iPoint, jPoint;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool jst_scheme = ((config->GetKind_Centered_AdjFlow() == CENTERED::JST) && (iMesh == MESH_0));
  bool grid_movement  = config->GetGrid_Movement();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge, normal, and neighbors---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));
    numerics->SetNeighbor(geometry->nodes->GetnNeighbor(iPoint), geometry->nodes->GetnNeighbor(jPoint));

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
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));
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
      Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock2Diag(jPoint, Jacobian_jj);
    }

  }
}


void CAdjEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  su2double Project_Grad_i, Project_Grad_j, *Limiter_i = nullptr,
  *Limiter_j = nullptr, *Psi_i = nullptr, *Psi_j = nullptr, *V_i, *V_j;
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;

  bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_AdjFlow() && (iMesh == MESH_0));
  bool limiter          = (config->GetKind_SlopeLimit_AdjFlow() != LIMITER::NONE);
  bool grid_movement    = config->GetGrid_Movement();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge and normal vectors ---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

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
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));
    }

    /*--- High order reconstruction using MUSCL strategy ---*/

    if (muscl) {

      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->nodes->GetCoord(jPoint, iDim) - geometry->nodes->GetCoord(iPoint, iDim));
        Vector_j[iDim] = 0.5*(geometry->nodes->GetCoord(iPoint, iDim) - geometry->nodes->GetCoord(jPoint, iDim));
      }

      /*--- Adjoint variables using gradient reconstruction and limiters ---*/

      const auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
      const auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

      if (limiter) {
        Limiter_i = nodes->GetLimiter(iPoint);
        Limiter_j = nodes->GetLimiter(jPoint);
      }

      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0; Project_Grad_j = 0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
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

      /* Check our reconstruction for exceeding bounds on the
       adjoint density. */

      su2double adj_limit = config->GetAdjointLimit();
      bool phi_bound_i = (fabs(Solution_i[0]) > adj_limit);
      bool phi_bound_j = (fabs(Solution_j[0]) > adj_limit);

      if (phi_bound_i) nodes->SetNon_Physical(iPoint, true);
      else             nodes->SetNon_Physical(iPoint, false);

      if (phi_bound_j) nodes->SetNon_Physical(jPoint, true);
      else             nodes->SetNon_Physical(jPoint, false);

      /* Lastly, check for existing first-order points still active
       from previous iterations. */

      if (nodes->GetNon_Physical(iPoint)) {
        counter_local++;
        for (iVar = 0; iVar < nVar; iVar++)
          Solution_i[iVar] = Psi_i[iVar];
      }
      if (nodes->GetNon_Physical(jPoint)) {
        counter_local++;
        for (iVar = 0; iVar < nVar; iVar++)
          Solution_j[iVar] = Psi_j[iVar];
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
      Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock2Diag(jPoint, Jacobian_jj);
    }

  }

  /*--- Warning message about non-physical reconstructions. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }

}

void CAdjEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  unsigned short iVar;
  unsigned long iPoint;
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  bool axisymmetric   = config->GetAxisymmetric();
  //  bool gravity        = (config->GetGravityForce() == YES);
  bool harmonic_balance  = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);

  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

  if (rotating_frame) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the adjoint variables ---*/
      numerics->SetAdjointVar(nodes->GetSolution(iPoint),
                              nodes->GetSolution(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the adjoint rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

    }
  }

  if (harmonic_balance) {

    su2double Volume, Source;

    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Get control volume ---*/
      Volume = geometry->nodes->GetVolume(iPoint);

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
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Set coordinate ---*/
      numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(iPoint));

      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);

      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

    }
  }

  //  if (gravity) {
  //
  //  }

}

void CAdjEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                      CConfig *config, unsigned short iMesh) {
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

    SharpEdge_Distance = (geometry->nodes->GetSharpEdge_Distance(iPoint) - config->GetAdjSharp_LimiterCoeff()*eps);

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

  SetResToZero();

  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->nodes->GetVolume(iPoint);
    Delta = solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) / Vol;

    Res_TruncError = nodes->GetResTruncError(iPoint);
    Residual = LinSysRes.GetBlock(iPoint);

    for (iVar = 0; iVar < nVar; iVar++) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      nodes->AddSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff);
      Residual_RMS[iVar] += Res*Res;
      AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
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

  SetResToZero();

  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->nodes->GetVolume(iPoint);
    Delta = solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) / Vol;

    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    local_Residual = LinSysRes.GetBlock(iPoint);

    for (iVar = 0; iVar < nVar; iVar++) {
      Res = local_Residual[iVar] + local_Res_TruncError[iVar];
      nodes->AddSolution(iPoint,iVar, -Res*Delta);
      Residual_RMS[iVar] += Res*Res;
      AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
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

  SetResToZero();

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/

    local_Res_TruncError = nodes->GetResTruncError(iPoint);

    /*--- Read the volume ---*/

    Vol = geometry->nodes->GetVolume(iPoint);

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
      Residual_RMS[iVar] += LinSysRes[total_index]*LinSysRes[total_index];
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
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
      nodes->AddSolution(iPoint,iVar, config->GetRelaxation_Factor_Adjoint()*LinSysSol[iPoint*nVar+iVar]);
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
  su2double *d = nullptr, *Normal = nullptr, *Psi = nullptr, *U = nullptr, Enthalpy, conspsi = 0.0, Mach_Inf,
  Area, *ConsPsi_Grad = nullptr,
  ConsPsi, d_press, grad_v, v_gradconspsi, UnitNormal[3], *GridVel = nullptr,
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
        if (geometry->nodes->GetDomain(iPoint)) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(iPoint);
          conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
          for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];

          nodes->SetAuxVar(iPoint,0,conspsi);

          /*--- Also load the auxiliary variable for first neighbors ---*/

          for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
            Neigh = geometry->nodes->GetPoint(iPoint, iNeigh);
            Psi = nodes->GetSolution(Neigh);
            U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(Neigh);
            Enthalpy = solver_container[FLOW_SOL]->GetNodes()->GetEnthalpy(Neigh);
            conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
            for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
            nodes->SetAuxVar(Neigh,0,conspsi);
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
        if (geometry->nodes->GetDomain(iPoint)) {

          d = nodes->GetForceProj_Vector(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = GeometryToolbox::Norm(nDim, Normal);

          const auto PrimVar_Grad = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
          ConsPsi_Grad = nodes->GetAuxVarGradient(iPoint)[0];
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
              GridVel = geometry->nodes->GetGridVel(iPoint);
              v_gradconspsi -= GridVel[iDim] * ConsPsi_Grad[iDim];
            }

          }

          /*--- Compute sensitivity for each surface point ---*/

          CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area * scale * factor;

          /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

          if (config->GetSens_Remove_Sharp()) {
            eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
            if ( geometry->nodes->GetSharpEdge_Distance(iPoint) < config->GetAdjSharp_LimiterCoeff()*eps )
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

        if (geometry->nodes->GetDomain(iPoint)) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          Area = GeometryToolbox::Norm(nDim, Normal);

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

        if (geometry->nodes->GetDomain(iPoint)) {
          Psi = nodes->GetSolution(iPoint);
          U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          r = U[0]; ru = U[1]; rv = U[2];
          if (nDim == 2) { rw = 0.0; rE = U[3]; }
          else { rw = U[3]; rE = U[4]; }
          p = Gamma_Minus_One*(rE-(ru*ru + rv*rv + rw*rw)/(2*r));

          Area = GeometryToolbox::Norm(nDim, Normal);
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

        if (geometry->nodes->GetDomain(iPoint)) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          p = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint);

          Mach_Inf   = config->GetMach();
          if (grid_movement) Mach_Inf = config->GetMach_Motion();

          d = nodes->GetForceProj_Vector(iPoint);
          Area = GeometryToolbox::Norm(nDim, Normal);
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

  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_BPress, &Total_Sens_BPress, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

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
        Coord_begin = geometry->nodes->GetCoord(iPoint);
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Coord_end = geometry->nodes->GetCoord(iPoint);
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
  su2double *d = nullptr, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitNormal;
  su2double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2;
  unsigned short iDim, iVar, jDim;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();

  UnitNormal = new su2double[nDim];
  Velocity = new su2double[nDim];
  Psi      = new su2double[nVar];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      /*--- Create a copy of the adjoint solution ---*/
      Psi_Aux = nodes->GetSolution(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];

      /*--- Flow solution ---*/
      U = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);

      /*--- Read the value of the objective function ---*/
      d = nodes->GetForceProj_Vector(iPoint);

      /*--- Normal vector computation ---*/
      Area = GeometryToolbox::Norm(nDim, Normal);
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
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
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
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
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
          su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel -= GridVel[iDim]*Normal[iDim];
          Jacobian_ii[0][0] -= ProjGridVel;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[iDim+1][iDim+1] -= ProjGridVel;
          Jacobian_ii[nVar-1][nVar-1] -= ProjGridVel;
        }

        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();

  Normal = new su2double[nDim];
  UnitNormal = new su2double[nDim];
  Velocity = new su2double[nDim];
  Psi_domain = new su2double[nVar];
  Psi_sym = new su2double[nVar];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

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
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
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
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
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
          su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel -= GridVel[iDim]*Normal[iDim];
          Jacobian_ii[0][0] -= ProjGridVel;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[iDim+1][iDim+1] -= ProjGridVel;
          Jacobian_ii[nVar-1][nVar-1] -= ProjGridVel;
        }

         /*--- Update jacobian ---*/

        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

      }

    }
  }

  delete [] Velocity;
  delete [] Psi_domain;
  delete [] Psi_sym;
  delete [] Normal;
  delete [] UnitNormal;

}

void CAdjEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint, Point_Normal;
  unsigned short iVar, iDim;
  su2double *Normal, *V_domain, *V_infty, *Psi_domain, *Psi_infty;

  bool implicit       = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_infty = new su2double[nVar];

  /*--- Loop over all the vertices ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- If the node belongs to the domain ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

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
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the upwind flux ---*/

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {

        /*--- Points in edge, coordinates and normal vector---*/
        su2double Coord_Reflected[3];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
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
          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

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
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {

        /*--- Index of the closest interior node ---*/

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Points in edge, coordinates and normal vector---*/

        su2double Coord_Reflected[3];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
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
          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_outlet = new su2double[nVar];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

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
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
      Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

      /*--- Viscous residual contribution (check again, Point_Normal was not being initialized before) ---*/

      if (config->GetViscous()) {

        /*--- Index of the closest interior node ---*/

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Points in edge, coordinates and normal vector---*/

        su2double Coord_Reflected[3];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
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
          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  INLET_TYPE Kind_Inlet = config->GetKind_Inlet();
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

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
      case INLET_TYPE::TOTAL_CONDITIONS:

        /*--- Adjoint solution at the inlet. Set to zero for now
             but should be replaced with derived expression for this type of
             inlet. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Psi_inlet[iVar] = 0.0;

        break;

        /*--- Mass flow has been specified at the inlet. ---*/
      case INLET_TYPE::MASS_FLOW:

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
          GridVel = geometry->nodes->GetGridVel(iPoint);
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
          bcn -= (1.0/Gamma_Minus_One)*ProjGridVel;
        }

        /*--- Impose value for PsiE based on hand-derived expression. ---*/
        Psi_inlet[nVar-1] = -phin*(1.0/bcn);

        break;

      default:
        SU2_MPI::Error("Unsupported INLET_TYPE.", CURRENT_FUNCTION);
        break;
      }


      /*--- Set the flow and adjoint states in the solver ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_inlet);

      /*--- Grid Movement ---*/

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {
        /*--- Index of the closest interior node ---*/

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Points in edge, coordinates and normal vector---*/

        su2double Coord_Reflected[3];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
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
          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
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

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Set the normal point ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      su2double Coord_Reflected[3];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      conv_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

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
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
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
          su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
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
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
            geometry->nodes->GetGridVel(iPoint));

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
          Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {

        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_outlet[nDim+5] = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        V_outlet[nDim+6] = solver_container[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint);

        /*--- Points in edge, coordinates and normal vector---*/
        visc_numerics->SetNormal(Normal);
        su2double Coord_Reflected[3];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

        /*--- Turbulent kinetic energy ---*/

        if (config->GetKind_Turb_Model() == TURB_MODEL::SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0), solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inflow = new su2double[nVar];

  /*--- Loop over all the vertices ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- If the node belong to the domain ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

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
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] Psi_domain;
  delete [] Psi_inflow;

}

void CAdjEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint, Point_Normal;
  su2double *Normal, *V_domain, *V_exhaust, *Psi_domain, *Psi_exhaust;
  unsigned short iVar, iDim;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_exhaust = new su2double[nVar];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

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
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

    }
  }

  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_exhaust;

}

void CAdjEulerSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  su2double *Normal, *V_domain, *V_inlet, *Psi_domain, *Psi_inlet;
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, GlobalIndex_inlet, GlobalIndex;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];

  /*--- Loop over all the vertices ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    GlobalIndex_inlet = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex = geometry->nodes->GetGlobalIndex(iPoint);

    /*--- If the node belong to the domain ---*/

    if ((geometry->nodes->GetDomain(iPoint)) && (GlobalIndex != GlobalIndex_inlet)) {

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
      su2double UnitNormal[3];

      Area = GeometryToolbox::Norm(nDim, Normal);

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
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar];
  Psi_outlet = new su2double[nVar];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    GlobalIndex_inlet = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex = geometry->nodes->GetGlobalIndex(iPoint);

    /*--- Check that the node belongs to the domain (i.e., not a halo node) and to discard the perimeter ---*/

    if ((geometry->nodes->GetDomain(iPoint)) && (GlobalIndex != GlobalIndex_inlet)) {

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
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);

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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool Grid_Movement = config->GetGrid_Movement();

  /*--- loop over points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Solution at time n-1, n and n+1 ---*/
    U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
    U_time_n   = nodes->GetSolution_time_n(iPoint);
    U_time_nP1 = nodes->GetSolution(iPoint);

    /*--- Volume at time n-1 and n ---*/
    if (Grid_Movement) {
      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_n = geometry->nodes->GetVolume_n(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);
    }
    else {
      Volume_nM1 = geometry->nodes->GetVolume(iPoint);
      Volume_n = geometry->nodes->GetVolume(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);
    }

    /*--- Time Step ---*/
    TimeStep = config->GetDelta_UnstTimeND();

    /*--- Compute Residual ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
        Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
      if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
        Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
                          +  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
    }

    /*--- Add Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual);

    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;

        if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
          Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
          Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
      }
      Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
    }
  }

}

void CAdjEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh;
  unsigned long index;
  su2double *Coord;
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
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, SU2_MPI::GetComm());
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
    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][ADJFLOW_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][ADJFLOW_SOL]->GetNodes()->GetSolution());
    solver[iMesh][ADJFLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][ADJFLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][ADJFLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;
  delete [] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

}

void CAdjEulerSolver::SetAuxVar_Surface_Gradient(CGeometry *geometry, const CConfig *config) {

  unsigned short iDim, jDim, iNeigh, iMarker;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint, iVertex;
  su2double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j;
  su2double **Smatrix, *Cvector;

  Smatrix = new su2double* [nDim];
  Cvector = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];


  /*--- Loop over boundary markers to select those for Euler or NS walls ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetSolid_Wall(iMarker)) {

      /*--- Loop over points on the surface (Least-Squares approximation) ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->nodes->GetDomain(iPoint)) {
          Coord_i = geometry->nodes->GetCoord(iPoint);
          AuxVar_i = nodes->GetAuxVar(iPoint);

          /*--- Inizialization of variables ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iDim] = 0.0;
          su2double r11 = 0.0, r12 = 0.0, r13 = 0.0, r22 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;

          for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
            jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
            Coord_j = geometry->nodes->GetCoord(jPoint);
            AuxVar_j = nodes->GetAuxVar(jPoint);

            su2double weight = 0;
            for (iDim = 0; iDim < nDim; iDim++)
              weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

            /*--- Sumations for entries of upper triangular matrix R ---*/
            r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
            r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
            r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
            if (nDim == 3) {
              r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
              r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
              r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
              r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
            }

            /*--- Entries of c:= transpose(A)*b ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/weight;
          }

          /*--- Entries of upper triangular matrix R ---*/
          r11 = sqrt(r11);
          r12 = r12/r11;
          r22 = sqrt(r22-r12*r12);
          if (nDim == 3) {
            r13 = r13/r11;
            r23 = r23_a/r22 - r23_b*r12/(r11*r22);
            r33 = sqrt(r33-r23*r23-r13*r13);
          }
          /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
          if (nDim == 2) {
            su2double detR2 = (r11*r22)*(r11*r22);
            Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
            Smatrix[0][1] = -r11*r12/detR2;
            Smatrix[1][0] = Smatrix[0][1];
            Smatrix[1][1] = r11*r11/detR2;
          }
          else {
            su2double detR2 = (r11*r22*r33)*(r11*r22*r33);
            su2double z11, z12, z13, z22, z23, z33; // aux vars
            z11 = r22*r33;
            z12 = -r12*r33;
            z13 = r12*r23-r13*r22;
            z22 = r11*r33;
            z23 = -r11*r23;
            z33 = r11*r22;
            Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
            Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
            Smatrix[0][2] = (z13*z33)/detR2;
            Smatrix[1][0] = Smatrix[0][1];
            Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
            Smatrix[1][2] = (z23*z33)/detR2;
            Smatrix[2][0] = Smatrix[0][2];
            Smatrix[2][1] = Smatrix[1][2];
            Smatrix[2][2] = (z33*z33)/detR2;
          }
          /*--- Computation of the gradient: S*c ---*/
          su2double product;
          for (iDim = 0; iDim < nDim; iDim++) {
            product = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              product += Smatrix[iDim][jDim]*Cvector[jDim];
            nodes->SetAuxVarGradient(iPoint, 0, iDim, product);
          }
        }
      } /*--- End of loop over surface points ---*/
    }
  }

  /*--- Memory deallocation ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Smatrix[iDim];
  delete [] Cvector;
  delete [] Smatrix;
}
