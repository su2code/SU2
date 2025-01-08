/*!
 * \file CRadP1Solver.cpp
 * \brief Main subroutines for solving P1 radiation problems.
 * \author Ruben Sanchez
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

#include "../../include/solvers/CRadP1Solver.hpp"
#include "../../include/variables/CRadP1Variable.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CRadP1Solver::CRadP1Solver() : CRadSolver() {

}

CRadP1Solver::CRadP1Solver(CGeometry* geometry, CConfig *config) : CRadSolver(geometry, config) {

  unsigned short iVar;
  unsigned short direct_diff = config->GetDirectDiff();
  bool multizone = config->GetMultizone_Problem();

  nDim =          geometry->GetnDim();
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
  nVar =          1;

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  Residual = new su2double[nVar];
  Solution = new su2double[nVar];

  Res_Visc = new su2double[nVar];

  /*--- Define some structures for locating max residuals ---*/

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Radiation() == EULER_IMPLICIT) {

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (P1 radiation equation)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

  }

  /*--- Solution and residual vectors ---*/

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Read farfield conditions from config ---*/
  Temperature_Inf = config->GetTemperature_FreeStreamND();

  /*--- Initialize the secondary values for direct derivative approxiations ---*/

  switch(direct_diff){
    case NO_DERIVATIVE: case D_DENSITY:
    case D_PRESSURE: case D_VISCOSITY:
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Not necessary here ---*/
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    default:
      break;
  }

  SetTemperature_Inf(Temperature_Inf);

  /*--- Initialize the BGS residuals. ---*/
  if (multizone){
    Residual_BGS.resize(nVar,1.0);
    Residual_Max_BGS.resize(nVar,1.0);
    Point_Max_BGS.resize(nVar,0);
    Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
  }

  /*--- Always instantiate and initialize the variable to a zero value. ---*/
  su2double init_val;
  switch(config->GetKind_P1_Init()){
    case P1_INIT::ZERO: init_val = 0.0; break;
    case P1_INIT::TEMPERATURE: init_val = 4.0*STEFAN_BOLTZMANN*pow(config->GetInc_Temperature_Init(),4.0); break;
    default: init_val = 0.0; break;
  }

  nodes = new CRadP1Variable(init_val, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Initialize the structure for volumetric heat source ---*/
  if(config->GetHeatSource()){
   SetVolumetricHeatSource(geometry, config);
  }

  SolverName = "RAD";
}

CRadP1Solver::~CRadP1Solver() {

  delete nodes;

}

void CRadP1Solver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;

  /*--- Initialize the residual vector ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    LinSysRes.SetBlock_Zero(iPoint);
  }

  /*--- Initialize the Jacobian matrix ---*/
  Jacobian.SetValZero();

  /*--- Compute the Solution gradients ---*/
  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config, true);
  }
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

}

void CRadP1Solver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  unsigned long iPoint;
  su2double Energy, Temperature;
  su2double SourceTerm, SourceTerm_Derivative;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Retrieve the radiative energy ---*/
    Energy = nodes->GetSolution(iPoint, 0);

    /*--- Retrieve temperature from the flow solver ---*/
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);

    /*--- Compute the divergence of the radiative flux ---*/
    SourceTerm = Absorption_Coeff*(Energy - 4.0*STEFAN_BOLTZMANN*pow(Temperature,4.0));

    /*--- Compute the derivative of the source term with respect to the temperature ---*/
    SourceTerm_Derivative =  - 16.0*Absorption_Coeff*STEFAN_BOLTZMANN*pow(Temperature,3.0);

    /*--- Store the source term and its derivative ---*/
    nodes->SetRadiative_SourceTerm(iPoint, 0, SourceTerm);
    nodes->SetRadiative_SourceTerm(iPoint, 1, SourceTerm_Derivative);

  }

}

void CRadP1Solver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[VISC_TERM];

  unsigned long iEdge, iPoint, jPoint;

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge ---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Points coordinates, and normal vector ---*/

    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                       geometry->nodes->GetCoord(jPoint));
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Radiation variables w/o reconstruction, and its gradients ---*/

    numerics->SetRadVar(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    numerics->SetRadVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(jPoint));

    /*--- Compute residual, and Jacobians ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Add and subtract residual, and update Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, Jacobian_i, Jacobian_j);

  }

}

void CRadP1Solver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                  CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), nullptr);

    /*--- Radiation variables w/o reconstruction ---*/

    numerics->SetRadVar(nodes->GetSolution(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Compute the source term ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

  }

}

void CRadP1Solver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {

  unsigned short iVar, jVar;
  unsigned long iVertex, iPoint;

  su2double Theta, Ib_w, Radiative_Energy;
  su2double *Normal, Area, Wall_Emissivity;
  su2double Radiative_Heat_Flux;
  su2double Twall;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall emissivity from config ---*/
  Wall_Emissivity = config->GetWall_Emissivity(Marker_Tag);

  /*--- Compute the constant for the wall theta ---*/
  Theta = Wall_Emissivity / (2.0*(2.0 - Wall_Emissivity));

    /*--- Retrieve the specified wall temperature ---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = GeometryToolbox::Norm(nDim, Normal);

      // Weak application of the boundary condition

      /*--- Initialize the viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Apply a weak boundary condition for the radiative transfer equation. ---*/

      /*--- Compute the blackbody intensity at the wall. ---*/
      Ib_w = 4.0*STEFAN_BOLTZMANN*pow(Twall,4.0);

      /*--- Compute the radiative heat flux. ---*/
      Radiative_Energy = nodes->GetSolution(iPoint, 0);
      Radiative_Heat_Flux = 1.0*Theta*(Ib_w - Radiative_Energy);

      /*--- Compute the Viscous contribution to the residual ---*/
      Res_Visc[0] = Radiative_Heat_Flux*Area;

      /*--- Apply to the residual vector ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Compute the Jacobian contribution. ---*/
      if (implicit) {
        Jacobian_i[0][0] = - Theta;
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
      }
    }
  }

}

void CRadP1Solver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iVar, jVar;
  unsigned long iVertex, iPoint;

  su2double Theta, Ib_w, Radiative_Energy;
  su2double *Normal, Area, Wall_Emissivity;
  su2double Radiative_Heat_Flux;
  su2double Twall;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall emissivity from config ---*/
  Wall_Emissivity = config->GetWall_Emissivity(Marker_Tag);

  /*--- Compute the constant for the wall theta ---*/
  Theta = Wall_Emissivity / (2.0*(2.0 - Wall_Emissivity));

  /*--- Retrieve the specified wall temperature ---*/
  Twall = GetTemperature_Inf();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = GeometryToolbox::Norm(nDim, Normal);

      // Weak application of the boundary condition

      /*--- Initialize the viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Apply a weak boundary condition for the radiative transfer equation. ---*/

      /*--- Compute the blackbody intensity at the wall. ---*/
      Ib_w = 4.0*STEFAN_BOLTZMANN*pow(Twall,4.0);

      /*--- Compute the radiative heat flux. ---*/
      Radiative_Energy = nodes->GetSolution(iPoint, 0);
      Radiative_Heat_Flux = 1.0*Theta*(Ib_w - Radiative_Energy);

      /*--- Compute the Viscous contribution to the residual ---*/
      Res_Visc[0] = Radiative_Heat_Flux*Area;

      /*--- Apply to the residual vector ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Compute the Jacobian contribution. ---*/
      if (implicit) {
        Jacobian_i[0][0] = - Theta;
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
      }
    }
  }

}

void CRadP1Solver::BC_Marshak(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, jVar;
  unsigned long iVertex, iPoint;

  su2double Theta, Ib_w, Temperature, Radiative_Energy;
  su2double *Normal, Area, Wall_Emissivity;
  su2double Radiative_Heat_Flux;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall emissivity from config ---*/
  Wall_Emissivity = config->GetWall_Emissivity(Marker_Tag);

  /*--- Compute the constant for the wall theta ---*/
  Theta = Wall_Emissivity / (2.0*(2.0 - Wall_Emissivity));

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = GeometryToolbox::Norm(nDim, Normal);

      // Weak application of the boundary condition

      /*--- Initialize the viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Apply a weak boundary condition for the radiative transfer equation. ---*/

      /*--- Retrieve temperature from the flow solver ---*/
      Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);

      /*--- Compute the blackbody intensity at the wall. ---*/
      Ib_w = 4.0*STEFAN_BOLTZMANN*pow(Temperature,4.0);

      /*--- Compute the radiative heat flux. ---*/
      Radiative_Energy = nodes->GetSolution(iPoint, 0);
      Radiative_Heat_Flux = Theta*(Ib_w - Radiative_Energy);

      /*--- Compute the Viscous contribution to the residual ---*/
      Res_Visc[0] = Radiative_Heat_Flux*Area;

      /*--- Apply to the residual vector ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Compute the Jacobian contribution. ---*/
      if (implicit) {
        Jacobian_i[0][0] = - Theta;
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
      }

    }
  }

}


void CRadP1Solver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Vol;
  su2double Delta;

  /*--- Set maximum residual to zero ---*/

  SetResToZero();

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the volume ---*/

    Vol = geometry->nodes->GetVolume(iPoint);

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    if (nodes->GetDelta_Time(iPoint) != 0.0) {
      Delta = Vol / nodes->GetDelta_Time(iPoint);
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index]);
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

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      nodes->AddSolution(iPoint, iVar, LinSysSol[iPoint*nVar+iVar]);
    }
  }

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CRadP1Solver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                               unsigned short iMesh, unsigned long Iteration) {

  unsigned short iMarker;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  su2double Area, Vol, Lambda;
  su2double Global_Delta_Time = 1E6, Local_Delta_Time = 0.0, K_v = 0.25;
  su2double CFL = config->GetCFL_Rad();
  su2double GammaP1 = 1.0 / (3.0*(Absorption_Coeff + Scattering_Coeff));
  const su2double* Normal;

  /*--- Compute spectral radius based on thermal conductivity ---*/

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetMax_Lambda_Visc(iPoint, 0.0);
  }

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Get the edge's normal vector to compute the edge's area ---*/
    Normal = geometry->edges->GetNormal(iEdge);
    Area = GeometryToolbox::Norm(nDim, Normal);

    /*--- Viscous contribution ---*/

    Lambda = GammaP1*Area*Area;
    if (geometry->nodes->GetDomain(iPoint)) nodes->AddMax_Lambda_Visc(iPoint, Lambda);
    if (geometry->nodes->GetDomain(jPoint)) nodes->AddMax_Lambda_Visc(jPoint, Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Viscous contribution ---*/

      Lambda = GammaP1*Area*Area;
      if (geometry->nodes->GetDomain(iPoint)) nodes->AddMax_Lambda_Visc(iPoint, Lambda);

    }
  }

  /*--- Each element uses their own speed, steady state simulation ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->nodes->GetVolume(iPoint);

    if (Vol != 0.0) {

      /*--- Time step setting method ---*/

       Local_Delta_Time = CFL*K_v*Vol*Vol/ nodes->GetMax_Lambda_Visc(iPoint);

      /*--- Min-Max-Logic ---*/

      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();

      nodes->SetDelta_Time(iPoint, Local_Delta_Time);
    }
    else {
      nodes->SetDelta_Time(iPoint, 0.0);
    }
  }

  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetComm_Level() == COMM_FULL) {

    su2double sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Allreduce(&sbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

    sbuf_time = Max_Delta_Time;
    SU2_MPI::Allreduce(&sbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  }

}
