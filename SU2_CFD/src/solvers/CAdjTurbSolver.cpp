/*!
 * \file CAdjTurbVariable.cpp
 * \brief Main subroutines for solving turbulent adjoint problems.
 * \author F. Palacios, A. Bueno, T. Economon
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


#include "../../include/solvers/CAdjTurbSolver.hpp"

CAdjTurbSolver::CAdjTurbSolver() : CSolver() {}

CAdjTurbSolver::CAdjTurbSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  adjoint = true;

  nDim = geometry->GetnDim();
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Dimension of the problem  ---*/
  nVar = config->GetnTurbVar();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar+1;

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  Residual   = new su2double [nVar];
  Residual_i = new su2double [nVar];
  Residual_j = new su2double [nVar];

  Solution   = new su2double [nVar];
  Solution_i = new su2double [nVar];
  Solution_j = new su2double [nVar];

  /*--- Define some auxiliar vector related with the geometry ---*/
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];

  /*--- Define some auxiliar vector related with the flow solution ---*/
  FlowSolution_i = new su2double [nDim+2]; FlowSolution_j = new su2double [nDim+2];

  /*--- Point to point Jacobians ---*/
  Jacobian_ii = new su2double* [nVar];
  Jacobian_ij = new su2double* [nVar];
  Jacobian_ji = new su2double* [nVar];
  Jacobian_jj = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_ii[iVar] = new su2double [nVar];
    Jacobian_ij[iVar] = new su2double [nVar];
    Jacobian_ji[iVar] = new su2double [nVar];
    Jacobian_jj[iVar] = new su2double [nVar];
  }

  /*--- Initialization of the structure of the whole Jacobian ---*/
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  Jacobian.SetValZero();
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Far-Field values and initizalization ---*/
  bool restart = config->GetRestart();

  nodes = new CAdjTurbVariable(0.0, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  if (restart && (iMesh == MESH_0)) {
    unsigned long index;
    su2double dull_val;
    string filename, AdjExt, text_line;
    ifstream restart_file;

    /*--- Restart the solution from file information ---*/
    string mesh_filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(mesh_filename);

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no file ---*/
    if (restart_file.fail()) {
      SU2_MPI::Error(string("There is no adjoint restart file ") + filename, CURRENT_FUNCTION);
    }

    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/

    getline (restart_file, text_line);

    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

      getline (restart_file, text_line);

      istringstream point_line(text_line);

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {

        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        nodes->SetSolution(iPoint_Local,0,Solution[0]);
        nodes->SetSolution_Old(iPoint_Local,0,Solution[0]);
        iPoint_Global_Local++;
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

    /*--- Close the restart file ---*/
    restart_file.close();
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

}

CAdjTurbSolver::~CAdjTurbSolver() {
  delete nodes;
}

void CAdjTurbSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;

  for (iVertex = 0; iVertex<geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      Solution[0] = 0.0;

      /*--- Set the solution values and zero the residual ---*/
      nodes->SetSolution_Old(iPoint,Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      Jacobian.DeleteValsRowi(iPoint);

    }
  }

}

void CAdjTurbSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;

  for (iVertex = 0; iVertex<geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      Solution[0] = 0.0;

      /*--- Set the solution values and zero the residual ---*/
      nodes->SetSolution_Old(iPoint,Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      Jacobian.DeleteValsRowi(iPoint);

    }
  }

}

void CAdjTurbSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Set Normal ---*/
    conv_numerics->SetNormal(geometry->vertex[val_marker][iVertex]->GetNormal());

    /*--- Set Conservative variables (for convection) ---*/
    su2double* U_i = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
    conv_numerics->SetConservative(U_i, nullptr);

    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    su2double* TurbPsi_i = nodes->GetSolution(iPoint);
    conv_numerics->SetTurbAdjointVar(TurbPsi_i, nullptr);

    /*--- Add Residuals and Jacobians ---*/
    conv_numerics->ComputeResidual(Residual, Jacobian_ii, nullptr, config);
    LinSysRes.AddBlock(iPoint, Residual);
    Jacobian.AddBlock2Diag(iPoint, Jacobian_ii);

  }

}

void CAdjTurbSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

  /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);

  }


    /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

  /*--- Gradient of the adjoint turbulent variables ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  /*--- Gradient of the turbulent variables ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solver_container[TURB_SOL]->SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solver_container[TURB_SOL]->SetSolution_Gradient_LS(geometry, config);

}

void CAdjTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  unsigned long iEdge, iPoint, jPoint;
  su2double *U_i, *U_j, *TurbPsi_i, *TurbPsi_j;

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge ---*/
    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
    U_j = solver_container[FLOW_SOL]->GetNodes()->GetSolution(jPoint);
    numerics->SetConservative(U_i, U_j);

    /*--- Set normal vectors and length ---*/
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    TurbPsi_i = nodes->GetSolution(iPoint);
    TurbPsi_j = nodes->GetSolution(jPoint);
    numerics->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);

    /*--- Gradient of turbulent variables w/o reconstruction ---*/
    const auto TurbVar_Grad_i = solver_container[TURB_SOL]->GetNodes()->GetGradient(iPoint);
    const auto TurbVar_Grad_j = solver_container[TURB_SOL]->GetNodes()->GetGradient(jPoint);
    numerics->SetScalarVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);

    /*--- Set normal vectors and length ---*/
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

    /*--- Add and Subtract Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);
    Jacobian.AddBlock2Diag(iPoint, Jacobian_ii);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
    Jacobian.AddBlock2Diag(jPoint, Jacobian_jj);

  }

}

void CAdjTurbSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                      CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[VISC_TERM];

  unsigned long iEdge, iPoint, jPoint;
  su2double *Coord_i, *Coord_j;

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge ---*/
    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Points coordinates, and set normal vectors and length ---*/
    Coord_i = geometry->nodes->GetCoord(iPoint);
    Coord_j = geometry->nodes->GetCoord(jPoint);
    numerics->SetCoord(Coord_i, Coord_j);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Conservative variables w/o reconstruction, turbulent variables w/o reconstruction,
     and turbulent adjoint variables w/o reconstruction ---*/
    numerics->SetConservative(solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint), solver_container[FLOW_SOL]->GetNodes()->GetSolution(jPoint));
    numerics->SetScalarVar(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint), solver_container[TURB_SOL]->GetNodes()->GetSolution(jPoint));
    numerics->SetTurbAdjointVar(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));

    /*--- Viscosity ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint),
                                  solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint));

    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    numerics->SetTurbAdjointGradient(nodes->GetGradient(iPoint), nodes->GetGradient(jPoint));

    /*--- Compute residual in a non-conservative way, and update ---*/
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

    /*--- Update adjoint viscous residual ---*/
    LinSysRes.AddBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);

    Jacobian.AddBlock2Diag(iPoint, Jacobian_ii);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
    Jacobian.AddBlock2Diag(jPoint, Jacobian_jj);

  }

}

void CAdjTurbSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];
  //CNumerics* second_numerics = numerics_container[SOURCE_SECOND_TERM];

  unsigned long iPoint;
  su2double *U_i, *TurbVar_i;
  su2double *TurbPsi_i; // Gradients

  /*--- Piecewise source term ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
    numerics->SetConservative(U_i, nullptr);

    /*--- Gradient of primitive variables w/o reconstruction ---*/
    auto GradPrimVar_i = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
    numerics->SetPrimVarGradient(GradPrimVar_i, nullptr);

    /*--- Laminar viscosity of the fluid ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint), 0.0);

    /*--- Turbulent variables w/o reconstruction ---*/
    TurbVar_i = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint);
    numerics->SetScalarVar(TurbVar_i, nullptr);

    /*--- Gradient of Turbulent Variables w/o reconstruction ---*/
    auto TurbVar_Grad_i = solver_container[TURB_SOL]->GetNodes()->GetGradient(iPoint);
    numerics->SetScalarVarGradient(TurbVar_Grad_i, nullptr);

    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    TurbPsi_i = nodes->GetSolution(iPoint);
    numerics->SetTurbAdjointVar(TurbPsi_i, nullptr);

    /*--- Gradient of Adjoint flow variables w/o reconstruction
     (for non-conservative terms depending on gradients of flow adjoint vars.) ---*/
    auto PsiVar_Grad_i = solver_container[ADJFLOW_SOL]->GetNodes()->GetGradient(iPoint);
    numerics->SetAdjointVarGradient(PsiVar_Grad_i, nullptr);

    /*--- Set volume and distances to the surface ---*/
    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));
    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Add and Subtract Residual ---*/
    numerics->ComputeResidual(Residual, Jacobian_ii, nullptr, config);
    LinSysRes.AddBlock(iPoint, Residual);
    Jacobian.AddBlock2Diag(iPoint, Jacobian_ii);

  }

//  /*--- Conservative Source Term ---*/
//  su2double **TurbVar_Grad_j;
//  unsigned long jPoint, iEdge;
//
//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//
//    /*--- Points in edge ---*/
//    iPoint = geometry->edges->GetNode(iEdge,0);
//    jPoint = geometry->edges->GetNode(iEdge,1);
//
//    /*--- Gradient of turbulent variables w/o reconstruction ---*/
//    TurbVar_Grad_i = solver_container[TURB_SOL]->GetNodes()->GetGradient(iPoint);
//    TurbVar_Grad_j = solver_container[TURB_SOL]->GetNodes()->GetGradient(jPoint);
//    second_numerics->SetScalarVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);
//
//    /*--- Turbulent adjoint variables w/o reconstruction ---*/
//    TurbPsi_i = nodes->GetSolution(iPoint);
//    TurbPsi_j = nodes->GetSolution(jPoint);
//    second_numerics->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);
//
//    /*--- Set normal vectors and length ---*/
//    second_numerics->SetNormal(geometry->edges->GetNormal(iEdge));
//
//    /*--- Add and Subtract Residual ---*/
//    second_numerics->ComputeResidual(Residual, Jacobian_ii, Jacobian_jj, config);
//    LinSysRes.AddBlock(iPoint, Residual);
//    LinSysRes.SubtractBlock(jPoint, Residual);
//    Jacobian.AddBlock2Diag(iPoint, Jacobian_ii);
//    Jacobian.AddBlock(iPoint, jPoint, Jacobian_jj);
//    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ii);
//    Jacobian.SubtractBlock2Diag(jPoint, Jacobian_jj);
//
//  }

}

void CAdjTurbSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol;

  /*--- Set maximum residual to zero ---*/

  SetResToZero();

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Read the volume ---*/

    Vol = geometry->nodes->GetVolume(iPoint);

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    Delta = Vol / (config->GetCFLRedCoeff_AdjTurb()*solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));

    Jacobian.AddVal2Diag(iPoint, Delta);

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = -LinSysRes[total_index];
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

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++)
      nodes->AddSolution(iPoint,iVar, LinSysSol[iPoint*nVar+iVar]);
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}
