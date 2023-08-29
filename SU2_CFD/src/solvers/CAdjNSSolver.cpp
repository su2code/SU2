/*!
 * \file CAdjNSSolver.cpp
 * \brief Main subroutines for solving Navier-Stokes adjoint problems.
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


#include "../../include/solvers/CAdjNSSolver.hpp"
#include "../../include/variables/CAdjNSVariable.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CAdjNSSolver::CAdjNSSolver() : CAdjEulerSolver() { }

CAdjNSSolver::CAdjNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CAdjEulerSolver() {
  unsigned long iPoint, iVertex;
  string text_line, mesh_filename;
  unsigned short iDim, iVar, iMarker;
  ifstream restart_file;
  string filename, AdjExt;

  adjoint = true;

  su2double RefArea    = config->GetRefArea();
  su2double RefDensity  = config->GetDensity_FreeStreamND();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double Mach_Motion     = config->GetMach_Motion();
  su2double *Normal = nullptr, myArea_Monitored;
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

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Define some auxiliary arrays related to the residual ---*/

  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]     = 0.0;
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

  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
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

CAdjNSSolver::~CAdjNSSolver() = default;


void CAdjNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh, unsigned long Iteration) {

  /*--- Use the flow solution to update the time step
   *    The time step depends on the characteristic velocity, which is the same
   *    for the adjoint and flow solutions, albeit in the opposite direction. ---*/
  solver_container[FLOW_SOL]->SetTime_Step(geometry, solver_container, config, iMesh, Iteration);

}

void CAdjNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  su2double SharpEdge_Distance;
  bool physical = true;

  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/

  bool implicit       = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool limiter        = (config->GetKind_SlopeLimit_AdjFlow() != LIMITER::NONE);
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

  /*--- Compute gradients adj for solution reconstruction and viscous term ---*/

  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetSolution_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
  }
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  /*--- Limiter computation (upwind reconstruction) ---*/

  if (limiter && !Output) SetSolution_Limiter(geometry, config);

  /*--- Compute gradients adj for viscous term coupling ---*/

  if ((config->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
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
    unsigned long MyErrorCounter = nonPhysicalPoints; nonPhysicalPoints = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &nonPhysicalPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(nonPhysicalPoints);
  }

}

void CAdjNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[VISC_TERM];

  unsigned long iPoint, jPoint, iEdge;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge, coordinates and normal vector---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(jPoint));
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

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
      Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.AddBlock2Diag(jPoint, Jacobian_jj);
    }

  }

}

void CAdjNSSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];
  CNumerics* second_numerics = numerics_container[SOURCE_SECOND_TERM];

  unsigned long iPoint, jPoint, iEdge;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();

  /*--- Loop over all the points, note that we are supposing that primitive and
   adjoint gradients have been computed previously ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Primitive variables w/o reconstruction, and its gradient ---*/

    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), nullptr);

    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint), nullptr);

    /*--- Gradient of adjoint variables ---*/

    numerics->SetAdjointVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/

    if ((config->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {

      /*--- Turbulent variables w/o reconstruction and its gradient ---*/

      numerics->SetScalarVar(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint), nullptr);

      numerics->SetScalarVarGradient(solver_container[TURB_SOL]->GetNodes()->GetGradient(iPoint), nullptr);

      /*--- Turbulent adjoint variables w/o reconstruction and its gradient ---*/

      numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->GetNodes()->GetSolution(iPoint), nullptr);

      numerics->SetTurbAdjointGradient(solver_container[ADJTURB_SOL]->GetNodes()->GetGradient(iPoint), nullptr);

      /*--- Set distance to the surface ---*/

      numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    }

    /*--- Compute residual ---*/

    numerics->ComputeResidual(Residual, config);

    /*--- Add to the residual ---*/

    LinSysRes.AddBlock(iPoint, Residual);

  }

  /*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/

  if ((config->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Points in edge, and normal vector ---*/

      iPoint = geometry->edges->GetNode(iEdge,0);
      jPoint = geometry->edges->GetNode(iEdge,1);
      second_numerics->SetNormal(geometry->edges->GetNormal(iEdge));

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

      second_numerics->SetScalarVar(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint),
                                solver_container[TURB_SOL]->GetNodes()->GetSolution(jPoint));

      /*--- Turbulent adjoint variables w/o reconstruction ---*/

      second_numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->GetNodes()->GetSolution(iPoint),
                                       solver_container[ADJTURB_SOL]->GetNodes()->GetSolution(jPoint));

      /*--- Set distance to the surface ---*/

      second_numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), geometry->nodes->GetWall_Distance(jPoint));

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
      second_numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the adjoint rotating frame source residual ---*/
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

    }
  }
}

void CAdjNSSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, jDim, iMarker, iPos, jPos;
  su2double *d = nullptr, div_phi, *Normal = nullptr, Area,
  normal_grad_psi5, normal_grad_T, sigma_partial, Laminar_Viscosity = 0.0, heat_flux_factor, temp_sens = 0.0,
  *Psi = nullptr, *U = nullptr, Enthalpy, gradPsi5_v, psi5_tau_partial, psi5_tau_grad_vel, source_v_1, Density,
  Pressure = 0.0, div_vel, val_turb_ke, vartheta, vartheta_partial, psi5_p_div_vel, Omega[3], rho_v[3] = {0.0,0.0,0.0},
  CrossProduct[3], r, ru, rv, rw, rE, p, T, dp_dr, dp_dru,
  dp_drv, dp_drw, dp_drE, dH_dr, dH_dru, dH_drv, dH_drw, dH_drE, H, D[3][3], Dd[3], Mach_Inf, eps, scale = 1.0,
  RefVel2, RefDensity, Mach2Vel, *Velocity_Inf, factor;

  auto *USens = new su2double[nVar];
  auto *UnitNormal = new su2double[nDim];
  auto *normal_grad_vel = new su2double[nDim];
  auto *tang_deriv_psi5 = new su2double[nDim];
  auto *tang_deriv_T = new su2double[nDim];
  auto **Sigma = new su2double* [nDim];

  for (iDim = 0; iDim < nDim; iDim++)
    Sigma[iDim] = new su2double [nDim];

  auto *normal_grad_gridvel = new su2double[nDim];
  auto *normal_grad_v_ux =new su2double[nDim];
  auto **Sigma_Psi5v = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Sigma_Psi5v[iDim] = new su2double [nDim];
  auto **tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    tau[iDim] = new su2double [nDim];
  auto *Velocity = new su2double[nDim];

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

        if (geometry->nodes->GetDomain(iPoint)) {

          const auto PsiVar_Grad = nodes->GetGradient(iPoint);
          const auto PrimVar_Grad = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);

          Laminar_Viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

          heat_flux_factor = Cp * Laminar_Viscosity / Prandtl_Lam;

          /*--- Compute face area and the unit normal to the surface ---*/

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = GeometryToolbox::Norm(nDim, Normal);
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
            // turb_ke is not considered in the stress tensor, see #797
            val_turb_ke = 0.0;
            CNumerics::ComputeStressTensor(nDim, tau, PrimVar_Grad+1, Laminar_Viscosity, Density, val_turb_ke);

            /*--- Form normal_grad_gridvel = \partial_n (u_omega) ---*/

            const auto GridVel_Grad = geometry->nodes->GetGridVel_Grad(iPoint);
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_grad_gridvel[iDim] = 0.0;
              for (jDim = 0; jDim < nDim; jDim++)
                normal_grad_gridvel[iDim] += GridVel_Grad[iDim][jDim]*UnitNormal[jDim];
            }

            /*--- Form normal_grad_v_ux = \partial_n (v - u_omega) ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              normal_grad_v_ux[iDim] = normal_grad_vel[iDim] - normal_grad_gridvel[iDim];

            /*--- Form Sigma_Psi5v ---*/

            div_vel = 0.0;
            for (iDim = 0 ; iDim < nDim; iDim++) {
              Velocity[iDim] = U[iDim+1]/Density;
              div_vel += PrimVar_Grad[iDim+1][iDim];
            }

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

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)) {

      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      Sens_BPress[iMarker] = 0.0;

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

  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

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

  su2double *d, l1psi, vartheta, Sigma_5, phi[3] = {};
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

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();

  auto *Psi = new su2double[nVar];
  auto **Tau = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Tau[iDim] = new su2double [nDim];
  auto *Velocity = new su2double[nDim];
  auto *Normal = new su2double[nDim];
  auto *Edge_Vector = new su2double[nDim];
  auto **GradPhi = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new su2double [nDim];
  auto *GradPsiE = new su2double [nDim];

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

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
        GridVel = geometry->nodes->GetGridVel(iPoint);
        for (iDim = 0; iDim < nDim; iDim++)
          phi[iDim] -= Psi[nDim+1]*GridVel[iDim];
      }

      /*--- Impose the value of the adjoint velocity as a strong boundary
       condition (Dirichlet). Fix the adjoint velocity and remove any addtional
       contribution to the residual at this node. ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        nodes->SetSolution_Old(iPoint,iDim+1, phi[iDim]);

      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes(iPoint, iDim+1) = 0.0;
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

        GridVel = geometry->nodes->GetGridVel(iPoint);
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

        const auto PsiVar_Grad = nodes->GetGradient(iPoint);
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

            Coord_i = geometry->nodes->GetCoord(iPoint);
            Coord_j = geometry->nodes->GetCoord(Point_Normal);
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

            Coord_i = geometry->nodes->GetCoord(iPoint);
            Coord_j = geometry->nodes->GetCoord(Point_Normal);
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
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);
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
  force_stress, Sigma_5, phi[3] = {0.0,0.0,0.0};
  su2double phis1, phis2, sq_vel, ProjVel, Enthalpy, *GridVel, phi_u, d_n;
  su2double Energy, ViscDens, XiDens, Density, SoundSpeed, Pressure, dPhiE_dn, Laminar_Viscosity, Eddy_Viscosity,
  Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  su2double kGTdotn=0.0, Area=0.0, Xi=0.0;

  auto *Psi = new su2double[nVar];
  auto **Tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Tau[iDim] = new su2double [nDim];
  auto *Velocity = new su2double[nDim];
  auto *Normal = new su2double[nDim];

  auto **GradPhi = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new su2double [nDim];
  auto *GradPsiE = new su2double [nDim];
  su2double *GradT;// = new su2double[nDim];
  su2double *GradP;
  su2double *GradDens;
  auto *dPoRho2 = new su2double[nDim];

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
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

    if (geometry->nodes->GetDomain(iPoint)) {

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
      Volume = geometry->nodes->GetVolume(iPoint);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Get the force projection vector (based on the objective function) ---*/
      d = nodes->GetForceProj_Vector(iPoint);

      /*--- Adjustments to strong boundary condition for dynamic meshes ---*/
      if ( grid_movement) {
        GridVel = geometry->nodes->GetGridVel(iPoint);
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
        LinSysRes(iPoint, iDim+1) = 0.0;
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

        Area = GeometryToolbox::Norm(nDim, Normal);

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
      LinSysRes(iPoint, nVar-1) = 0.0;
      nodes->SetEnergy_ResTruncError_Zero(iPoint);
      nodes->SetSolution_Old(iPoint,nDim+1, q);
      if (implicit) {
        iVar = nDim+1;
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

      /*--- Additional contributions to adjoint density (weak imposition) ---*/

      /*--- Acquire gradient information ---*/
      auto PsiVar_Grad = nodes->GetGradient(iPoint);
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
        GridVel = geometry->nodes->GetGridVel(iPoint);

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
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_ii);
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
