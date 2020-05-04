/*!
 * \file CDiscAdjFEASolver.cpp
 * \brief Main subroutines for solving adjoint FEM elasticity problems.
 * \author R. Sanchez
 * \version 7.0.4 "Blackbird"
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


#include "../../include/solvers/CDiscAdjFEASolver.hpp"
#include "../../include/variables/CDiscAdjFEAVariable.hpp"

CDiscAdjFEASolver::CDiscAdjFEASolver(void) : CSolver() { }

CDiscAdjFEASolver::CDiscAdjFEASolver(CGeometry *geometry, CConfig *config)  : CSolver() { }

CDiscAdjFEASolver::CDiscAdjFEASolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver,
                                     unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  adjoint = true;

  unsigned short iVar, iMarker;

  unsigned long iPoint;
  string text_line, mesh_filename;
  string filename, AdjExt;

  bool dynamic = (config->GetTime_Domain());

  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max = new unsigned long[nVar]();
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim]();
  }

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()) {

    Residual_BGS      = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]      = 1.0;
    Residual_Max_BGS  = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS = new unsigned long[nVar]();
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim]();
    }

  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 1e-16;

  if (dynamic) {
    Solution_Vel    = new su2double[nVar];
    Solution_Accel  = new su2double[nVar];

    for (iVar = 0; iVar < nVar; iVar++) Solution_Vel[iVar]      = 1e-16;
    for (iVar = 0; iVar < nVar; iVar++) Solution_Accel[iVar]    = 1e-16;
  }

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]]();
  }

  Sens_E  = new su2double[nMarker]();
  Sens_Nu = new su2double[nMarker]();
  Sens_nL = new su2double[nMarker]();

  nodes = new CDiscAdjFEABoundVariable(Solution, Solution_Accel, Solution_Vel, nPoint, nDim, nVar, dynamic, config);
  SetBaseClassPointerToNodes();

  /*--- Set which points are vertices and allocate boundary data. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
      if (iVertex >= 0) {
        nodes->Set_isVertex(iPoint,true);
        break;
      }
    }
  nodes->AllocateBoundaryVariables(config);


  /*--- Store the direct solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    nodes->SetSolution_Direct(iPoint, direct_solver->GetNodes()->GetSolution(iPoint));
  }

  if (dynamic){
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      nodes->SetSolution_Accel_Direct(iPoint, direct_solver->GetNodes()->GetSolution_Accel(iPoint));
    }

    for (iPoint = 0; iPoint < nPoint; iPoint++){
      nodes->SetSolution_Vel_Direct(iPoint, direct_solver->GetNodes()->GetSolution_Vel(iPoint));
    }
  }

  /*--- Initialize vector structures for multiple material definition ---*/

  nMPROP = config->GetnElasticityMod();

  /*--- For a material to be fully defined, we need to have the same number for all three parameters ---*/
  bool checkDef = ((config->GetnElasticityMod() == config->GetnPoissonRatio()) &&
                   (config->GetnElasticityMod() == config->GetnMaterialDensity()) &&
                   (config->GetnMaterialDensity() == config->GetnPoissonRatio()));

  if (!checkDef){
    SU2_MPI::Error("WARNING: For a material to be fully defined, E, Nu and Rho need to have the same dimensions.", CURRENT_FUNCTION);
  }

  E_i           = new su2double[nMPROP]();
  Local_Sens_E  = new su2double[nMPROP]();
  Global_Sens_E = new su2double[nMPROP]();
  Total_Sens_E  = new su2double[nMPROP]();
  AD_Idx_E_i    = new int[nMPROP]();

  Nu_i           = new su2double[nMPROP]();
  Local_Sens_Nu  = new su2double[nMPROP]();
  Global_Sens_Nu = new su2double[nMPROP]();
  Total_Sens_Nu  = new su2double[nMPROP]();
  AD_Idx_Nu_i    = new int[nMPROP]();

  Rho_i           = new su2double[nMPROP](); // For inertial effects
  Local_Sens_Rho  = new su2double[nMPROP]();
  Global_Sens_Rho = new su2double[nMPROP]();
  Total_Sens_Rho  = new su2double[nMPROP]();
  AD_Idx_Rho_i    = new int[nMPROP]();

  Rho_DL_i           = new su2double[nMPROP](); // For dead loads
  Local_Sens_Rho_DL  = new su2double[nMPROP]();
  Global_Sens_Rho_DL = new su2double[nMPROP]();
  Total_Sens_Rho_DL  = new su2double[nMPROP]();
  AD_Idx_Rho_DL_i    = new int[nMPROP]();

  /*--- Initialize vector structures for multiple electric regions ---*/

  de_effects = config->GetDE_Effects();

  if (de_effects) {
    nEField = config->GetnElectric_Field();

    EField             = new su2double[nEField]();
    Local_Sens_EField  = new su2double[nEField]();
    Global_Sens_EField = new su2double[nEField]();
    Total_Sens_EField  = new su2double[nEField]();
    AD_Idx_EField      = new int[nEField]();
  }

  /*--- Initialize vector structures for structural-based design variables ---*/

  switch (config->GetDV_FEA()) {
    case YOUNG_MODULUS:
    case POISSON_RATIO:
    case DENSITY_VAL:
    case DEAD_WEIGHT:
    case ELECTRIC_FIELD:
      fea_dv = true;
      break;
    default:
      fea_dv = false;
      break;
  }

  if (fea_dv) {
    ReadDV(config);
    Local_Sens_DV  = new su2double[nDV]();
    Global_Sens_DV = new su2double[nDV]();
    Total_Sens_DV  = new su2double[nDV]();
    AD_Idx_DV_Val  = new int[nDV]();
  }

}

CDiscAdjFEASolver::~CDiscAdjFEASolver(void){

  unsigned short iMarker;

  if (CSensitivity != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

  delete [] E_i;
  delete [] Nu_i;
  delete [] Rho_i;
  delete [] Rho_DL_i;

  delete [] AD_Idx_E_i;
  delete [] AD_Idx_Nu_i;
  delete [] AD_Idx_Rho_i;
  delete [] AD_Idx_Rho_DL_i;

  delete [] Local_Sens_E;
  delete [] Local_Sens_Nu;
  delete [] Local_Sens_Rho;
  delete [] Local_Sens_Rho_DL;

  delete [] Global_Sens_E;
  delete [] Global_Sens_Nu;
  delete [] Global_Sens_Rho;
  delete [] Global_Sens_Rho_DL;

  delete [] Total_Sens_E;
  delete [] Total_Sens_Nu;
  delete [] Total_Sens_Rho;
  delete [] Total_Sens_Rho_DL;

  delete [] normalLoads;
  delete [] Sens_E;
  delete [] Sens_Nu;
  delete [] Sens_nL;

  delete [] EField;
  delete [] Local_Sens_EField;
  delete [] Global_Sens_EField;
  delete [] Total_Sens_EField;
  delete [] AD_Idx_EField;

  delete [] DV_Val;
  delete [] Local_Sens_DV;
  delete [] Global_Sens_DV;
  delete [] Total_Sens_DV;
  delete [] AD_Idx_DV_Val;

  delete [] Solution_Vel;
  delete [] Solution_Accel;

  delete nodes;
}

void CDiscAdjFEASolver::SetRecording(CGeometry* geometry, CConfig *config){


  bool dynamic (config->GetTime_Domain());

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint, nodes->GetSolution_Direct(iPoint));
  }

  if (dynamic){
    /*--- Reset the solution to the initial (converged) solution ---*/

    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->GetNodes()->SetSolution_Accel(iPoint, nodes->GetSolution_Accel_Direct(iPoint));
    }

    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->GetNodes()->SetSolution_Vel(iPoint, nodes->GetSolution_Vel_Direct(iPoint));
    }

    /*--- Reset the input for time n ---*/

    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n(iPoint)[iVar]);
      }
    }
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_Accel_time_n(iPoint)[iVar]);
      }
    }
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_Vel_time_n(iPoint)[iVar]);
      }
    }

  }

  /*--- Set the Jacobian to zero since this is not done inside the meanflow iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjFEASolver::RegisterSolution(CGeometry *geometry, CConfig *config){

  bool input = true;
  bool dynamic = config->GetTime_Domain();
  bool push_index = !config->GetMultizone_Problem();

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  direct_solver->GetNodes()->RegisterSolution(input, push_index);

  if (dynamic) {

    /*--- Register acceleration (u'') and velocity (u') at time step n ---*/

    direct_solver->GetNodes()->RegisterSolution_Accel(input);
    direct_solver->GetNodes()->RegisterSolution_Vel(input);

    /*--- Register solution (u), acceleration (u'') and velocity (u') at time step n-1 ---*/

    direct_solver->GetNodes()->Register_femSolution_time_n();
    direct_solver->GetNodes()->RegisterSolution_Accel_time_n();
    direct_solver->GetNodes()->RegisterSolution_Vel_time_n();

  }

}

void CDiscAdjFEASolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){

  /*--- Register element-based values as input ---*/

  unsigned short iVar;

  if (KindDirect_Solver == RUNTIME_FEA_SYS) {

    bool pseudo_static = config->GetPseudoStatic();

    for (iVar = 0; iVar < nMPROP; iVar++) {
      E_i[iVar]      = config->GetElasticyMod(iVar);
      Nu_i[iVar]     = config->GetPoissonRatio(iVar);
      Rho_i[iVar]    = pseudo_static? 0.0 : config->GetMaterialDensity(iVar);
      Rho_DL_i[iVar] = config->GetMaterialDensity(iVar);
    }

    /*--- Read the values of the electric field ---*/
    if (de_effects) {
      for (iVar = 0; iVar < nEField; iVar++)
        EField[iVar] = config->Get_Electric_Field_Mod(iVar);
    }

    /*--- Reset index, otherwise messes up other derivatives ---*/
    if (fea_dv) {
      for (iVar = 0; iVar < nDV; iVar++) AD::ResetInput(DV_Val[iVar]);
    }

    if (!reset) {
      bool local_index = config->GetMultizone_Problem();
      bool push_index = !local_index;

      for (iVar = 0; iVar < nMPROP; iVar++) {
        AD::RegisterInput(E_i[iVar], push_index);
        AD::RegisterInput(Nu_i[iVar], push_index);
        AD::RegisterInput(Rho_i[iVar], push_index);
        AD::RegisterInput(Rho_DL_i[iVar], push_index);
      }

      if(de_effects){
        for (iVar = 0; iVar < nEField; iVar++)
          AD::RegisterInput(EField[iVar], push_index);
      }

      if(fea_dv){
        for (iVar = 0; iVar < nDV; iVar++)
          AD::RegisterInput(DV_Val[iVar], push_index);
      }

      /*--- Explicitly store the tape indices for when we extract the derivatives ---*/
      if (local_index) {
        for (iVar = 0; iVar < nMPROP; iVar++) {
          AD::SetIndex(AD_Idx_E_i[iVar], E_i[iVar]);
          AD::SetIndex(AD_Idx_Nu_i[iVar], Nu_i[iVar]);
          AD::SetIndex(AD_Idx_Rho_i[iVar], Rho_i[iVar]);
          AD::SetIndex(AD_Idx_Rho_DL_i[iVar], Rho_DL_i[iVar]);
        }

        if (de_effects) {
          for (iVar = 0; iVar < nEField; iVar++)
            AD::SetIndex(AD_Idx_EField[iVar], EField[iVar]);
        }

        if (fea_dv) {
          for (iVar = 0; iVar < nDV; iVar++)
            AD::SetIndex(AD_Idx_DV_Val[iVar], DV_Val[iVar]);
        }
      }

      /*--- Register the flow tractions ---*/
      if (config->GetnMarker_Fluid_Load() > 0)
        direct_solver->GetNodes()->RegisterFlowTraction();
    }

  }

    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/
}

void CDiscAdjFEASolver::RegisterOutput(CGeometry *geometry, CConfig *config){

  bool input = false;
  bool dynamic = config->GetTime_Domain();
  bool push_index = !config->GetMultizone_Problem();

  /*--- Register variables as output of the solver iteration ---*/

  direct_solver->GetNodes()->RegisterSolution(input, push_index);

  if (dynamic) {
    /*--- Register acceleration (u'') and velocity (u') at time step n ---*/
    direct_solver->GetNodes()->RegisterSolution_Accel(input);
    direct_solver->GetNodes()->RegisterSolution_Vel(input);
  }

}

void CDiscAdjFEASolver::RegisterObj_Func(CConfig *config){

  /*--- Here we can add new (scalar) objective functions ---*/

  switch (config->GetKind_ObjFunc()){
  case REFERENCE_GEOMETRY:
      ObjFunc_Value = direct_solver->GetTotal_OFRefGeom();
      break;
  case REFERENCE_NODE:
      ObjFunc_Value = direct_solver->GetTotal_OFRefNode();
      break;
  case VOLUME_FRACTION:
  case TOPOL_DISCRETENESS:
      ObjFunc_Value = direct_solver->GetTotal_OFVolFrac();
      break;
  case TOPOL_COMPLIANCE:
      ObjFunc_Value = direct_solver->GetTotal_OFCompliance();
      break;
  default:
      ObjFunc_Value = 0.0;  // If the objective function is computed in a different physical problem
      break;
 /*--- Template for new objective functions where TemplateObjFunction()
  *  is the routine that returns the obj. function value. The computation
  * must be done while the tape is active, i.e. between AD::StartRecording() and
  * AD::StopRecording() in DiscAdjMeanFlowIteration::Iterate(). The best place is somewhere
  * inside MeanFlowIteration::Iterate().
  *
  * case TEMPLATE_OBJECTIVE:
  *    ObjFunc_Value = TemplateObjFunction();
  *    break;
  * ---*/
  }
  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc_Value);
  }
}


void CDiscAdjFEASolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config){

  bool dynamic = (config->GetTime_Domain());
  unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
  unsigned long TimeIter = config->GetTimeIter();
  su2double seeding = 1.0;

  if (dynamic){
    if (TimeIter < IterAvg_Obj){
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else{
      seeding = 0.0;
    }
  }

  if (rank == MASTER_NODE){
    SU2_TYPE::SetDerivative(ObjFunc_Value, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
  }
}

void CDiscAdjFEASolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  bool dynamic = config->GetTime_Domain();
  bool multizone = config->GetMultizone_Problem();

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
    SetRes_RMS(iVar,0.0);
    SetRes_Max(iVar,0.0,0);
  }

  /*--- Set the old solution, for multi-zone problems this is done after computing the
   *    residuals, otherwise the per-zone-residuals do not make sense, as on entry Solution
   *    contains contributions from other zones but on extraction it does not. ---*/

  if(!multizone) nodes->Set_OldSolution();

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution ---*/

    if(config->GetMultizone_Problem()) {
      direct_solver->GetNodes()->GetAdjointSolution_LocalIndex(iPoint,Solution);
    }
    else {
      direct_solver->GetNodes()->GetAdjointSolution(iPoint,Solution);
    }

    /*--- Store the adjoint solution ---*/

    nodes->SetSolution(iPoint,Solution);

  }

  /*--- Solution for acceleration (u'') and velocity (u') at time n ---*/

  if (dynamic){

    /*--- FIRST: The acceleration solution ---*/

    /*--- Set the old acceleration solution ---*/
    nodes->Set_OldSolution_Accel();

    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint acceleration solution u'' ---*/

      direct_solver->GetNodes()->GetAdjointSolution_Accel(iPoint,Solution_Accel);

      /*--- Store the adjoint acceleration solution u'' ---*/

      nodes->SetSolution_Accel(iPoint,Solution_Accel);

    }

    /*--- NEXT: The velocity solution ---*/

    /*--- Set the old velocity solution ---*/
    nodes->Set_OldSolution_Vel();

    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint velocity solution u'' ---*/

      direct_solver->GetNodes()->GetAdjointSolution_Vel(iPoint,Solution_Vel);

      /*--- Store the adjoint velocity solution u'' ---*/

      nodes->SetSolution_Vel(iPoint,Solution_Vel);

    }

    /*--- NOW: The solution at time n ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint solution at time n ---*/

      direct_solver->GetNodes()->GetAdjointSolution_time_n(iPoint,Solution);

      /*--- Store the adjoint solution at time n ---*/

      nodes->Set_Solution_time_n(iPoint,Solution);
    }

    /*--- The acceleration solution at time n... ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint acceleration solution u'' at time n ---*/

      direct_solver->GetNodes()->GetAdjointSolution_Accel_time_n(iPoint,Solution_Accel);

      /*--- Store the adjoint acceleration solution u'' at time n---*/

      nodes->SetSolution_Accel_time_n(iPoint,Solution_Accel);

    }

    /*--- ... and the velocity solution at time n ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint velocity solution u' at time n ---*/

      direct_solver->GetNodes()->GetAdjointSolution_Vel_time_n(iPoint,Solution_Vel);

      /*--- Store the adjoint velocity solution u' at time n ---*/

      nodes->SetSolution_Vel_time_n(iPoint,Solution_Vel);

    }

  }

  /*--- TODO: Need to set the MPI solution in the previous TS ---*/

  /*--- Set the residuals ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      residual = nodes->GetSolution(iPoint, iVar) - nodes->GetSolution_Old(iPoint, iVar);

      AddRes_RMS(iVar,residual*residual);
      AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
    }
    if (dynamic){
      for (iVar = 0; iVar < nVar; iVar++){
        residual = nodes->GetSolution_Accel(iPoint, iVar) - nodes->GetSolution_Old_Accel(iPoint, iVar);

        AddRes_RMS(iVar,residual*residual);
        AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
      for (iVar = 0; iVar < nVar; iVar++){
        residual = nodes->GetSolution_Vel(iPoint, iVar) - nodes->GetSolution_Old_Vel(iPoint, iVar);

        AddRes_RMS(iVar,residual*residual);
        AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
    }
  }

  SetResidual_RMS(geometry, config);

}

void CDiscAdjFEASolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  bool local_index = config->GetMultizone_Problem();

  /*--- Extract the adjoint values of the farfield values ---*/

  if (KindDirect_Solver == RUNTIME_FEA_SYS){

    if (local_index) {
      for (iVar = 0; iVar < nMPROP; iVar++) {
        Local_Sens_E[iVar] = AD::GetDerivative(AD_Idx_E_i[iVar]);
        Local_Sens_Nu[iVar] = AD::GetDerivative(AD_Idx_Nu_i[iVar]);
        Local_Sens_Rho[iVar] = AD::GetDerivative(AD_Idx_Rho_i[iVar]);
        Local_Sens_Rho_DL[iVar] = AD::GetDerivative(AD_Idx_Rho_DL_i[iVar]);
      }
    }
    else {
      for (iVar = 0; iVar < nMPROP; iVar++) {
        Local_Sens_E[iVar] = SU2_TYPE::GetDerivative(E_i[iVar]);
        Local_Sens_Nu[iVar] = SU2_TYPE::GetDerivative(Nu_i[iVar]);
        Local_Sens_Rho[iVar] = SU2_TYPE::GetDerivative(Rho_i[iVar]);
        Local_Sens_Rho_DL[iVar] = SU2_TYPE::GetDerivative(Rho_DL_i[iVar]);
      }
    }

    SU2_MPI::Allreduce(Local_Sens_E, Global_Sens_E,  nMPROP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Local_Sens_Nu, Global_Sens_Nu, nMPROP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Local_Sens_Rho, Global_Sens_Rho, nMPROP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Local_Sens_Rho_DL, Global_Sens_Rho_DL, nMPROP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /*--- Extract the adjoint values of the electric field in the case that it is a parameter of the problem. ---*/

    if (de_effects) {
      for (iVar = 0; iVar < nEField; iVar++) {
        if (local_index) Local_Sens_EField[iVar] = AD::GetDerivative(AD_Idx_EField[iVar]);
        else             Local_Sens_EField[iVar] = SU2_TYPE::GetDerivative(EField[iVar]);
      }
      SU2_MPI::Allreduce(Local_Sens_EField, Global_Sens_EField, nEField, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    if (fea_dv) {
      for (iVar = 0; iVar < nDV; iVar++) {
        if (local_index) Local_Sens_DV[iVar] = AD::GetDerivative(AD_Idx_DV_Val[iVar]);
        else             Local_Sens_DV[iVar] = SU2_TYPE::GetDerivative(DV_Val[iVar]);
      }
      SU2_MPI::Allreduce(Local_Sens_DV, Global_Sens_DV, nDV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    /*--- Extract the flow traction sensitivities ---*/

    if (config->GetnMarker_Fluid_Load() > 0){
      su2double val_sens;
      for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iDim = 0; iDim < nDim; iDim++){
          val_sens = direct_solver->GetNodes()->ExtractFlowTraction_Sensitivity(iPoint,iDim);
          nodes->SetFlowTractionSensitivity(iPoint, iDim, val_sens);
        }
      }
    }

  }

}

void CDiscAdjFEASolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config){

  bool dynamic = (config->GetTime_Domain());
  bool deform_mesh = (config->GetnMarker_Deform_Mesh() > 0);

  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      Solution[iVar] = nodes->GetSolution(iPoint,iVar);
    }
    if(deform_mesh){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] += nodes->GetSourceTerm_DispAdjoint(iPoint,iVar);
      }
    }
    if (dynamic){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution_Accel[iVar] = nodes->GetSolution_Accel(iPoint,iVar);
      }
      for (iVar = 0; iVar < nVar; iVar++){
        Solution_Vel[iVar] = nodes->GetSolution_Vel(iPoint,iVar);
      }
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] += nodes->GetDynamic_Derivative_n(iPoint,iVar);
      }
      for (iVar = 0; iVar < nVar; iVar++){
        Solution_Accel[iVar] += nodes->GetDynamic_Derivative_Accel_n(iPoint,iVar);
      }
      for (iVar = 0; iVar < nVar; iVar++){
        Solution_Vel[iVar] += nodes->GetDynamic_Derivative_Vel_n(iPoint,iVar);
      }
    }
    direct_solver->GetNodes()->SetAdjointSolution(iPoint,Solution);
    if (dynamic){
      direct_solver->GetNodes()->SetAdjointSolution_Accel(iPoint,Solution_Accel);
      direct_solver->GetNodes()->SetAdjointSolution_Vel(iPoint,Solution_Vel);
    }

  }

}

void CDiscAdjFEASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output){

  bool dynamic = (config_container->GetTime_Domain());
  unsigned long iPoint;
  unsigned short iVar;

  if (dynamic){
    for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++){
      for (iVar=0; iVar < nVar; iVar++){
        nodes->SetDynamic_Derivative_n(iPoint, iVar, nodes->GetSolution_time_n(iPoint, iVar));
      }
      for (iVar=0; iVar < nVar; iVar++){
        nodes->SetDynamic_Derivative_Accel_n(iPoint, iVar, nodes->GetSolution_Accel_time_n(iPoint, iVar));
      }
      for (iVar=0; iVar < nVar; iVar++){
        nodes->SetDynamic_Derivative_Vel_n(iPoint, iVar, nodes->GetSolution_Vel_time_n(iPoint, iVar));
      }
    }
  }

}

void CDiscAdjFEASolver::SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config){

  unsigned short iVar;

  for (iVar = 0; iVar < nMPROP; iVar++){
    Total_Sens_E[iVar]        += Global_Sens_E[iVar];
    Total_Sens_Nu[iVar]       += Global_Sens_Nu[iVar];
    Total_Sens_Rho[iVar]      += Global_Sens_Rho[iVar];
    Total_Sens_Rho_DL[iVar]   += Global_Sens_Rho_DL[iVar];
  }

  if (de_effects){
    for (iVar = 0; iVar < nEField; iVar++)
      Total_Sens_EField[iVar]+= Global_Sens_EField[iVar];
  }

  if (fea_dv){
    for (iVar = 0; iVar < nDV; iVar++)
      Total_Sens_DV[iVar] += Global_Sens_DV[iVar];
  }

  /*--- Extract the topology optimization density sensitivities ---*/

  direct_solver->ExtractAdjoint_Variables(geometry, config);

  /*--- Extract the geometric sensitivities ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    su2double *Coord = geometry->node[iPoint]->GetCoord();

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {

      su2double Sensitivity;

      if(config->GetMultizone_Problem()) {
        Sensitivity = geometry->node[iPoint]->GetAdjointSolution(iDim);
      }
      else {
        Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);
        /*--- Set the index manually to zero. ---*/
        AD::ResetInput(Coord[iDim]);
      }

      nodes->SetSensitivity(iPoint, iDim, Sensitivity);
    }
  }
  SetSurface_Sensitivity(geometry, config);
}

void CDiscAdjFEASolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config){

}

void CDiscAdjFEASolver::ReadDV(CConfig *config) {

  unsigned long index;

  string filename;
  ifstream properties_file;

  /*--- Choose the filename of the design variable ---*/

  string input_name;

  switch (config->GetDV_FEA()) {
    case YOUNG_MODULUS:
      input_name = "dv_young.opt";
      break;
    case POISSON_RATIO:
      input_name = "dv_poisson.opt";
      break;
    case DENSITY_VAL:
    case DEAD_WEIGHT:
      input_name = "dv_density.opt";
      break;
    case ELECTRIC_FIELD:
      input_name = "dv_efield.opt";
      break;
    default:
      input_name = "dv.opt";
      break;
  }

  filename = input_name;

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {

    if (rank == MASTER_NODE)
      cout << "There is no design variable file." << endl;

    nDV   = 1;
    DV_Val = new su2double[nDV];
    for (unsigned short iDV = 0; iDV < nDV; iDV++)
      DV_Val[iDV] = 1.0;

  }
  else{

    string text_line;

     /*--- First pass: determine number of design variables ---*/

    unsigned short iDV = 0;

    /*--- Skip the first line: it is the header ---*/

    getline (properties_file, text_line);

    while (getline (properties_file, text_line)) iDV++;

    /*--- Close the restart file ---*/

    properties_file.close();

    nDV = iDV;
    DV_Val = new su2double[nDV];

    /*--- Reopen the file (TODO: improve this) ---*/

    properties_file.open(filename.data(), ios::in);

    /*--- Skip the first line: it is the header ---*/

    getline (properties_file, text_line);

    iDV = 0;
    while (getline (properties_file, text_line)) {

      istringstream point_line(text_line);

      point_line >> index >> DV_Val[iDV];

      iDV++;

    }

    /*--- Close the restart file ---*/

    properties_file.close();

  }

}

void CDiscAdjFEASolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  unsigned short iVar;
  unsigned long index, counter;
  string restart_filename, filename;

  /*--- Restart the solution from file information ---*/

  filename = config->GetSolution_AdjFileName();
  restart_filename = config->GetObjFunc_Extension(filename);
  restart_filename = config->GetFilename(restart_filename, "", val_iter);

  /*--- Read and store the restart metadata. ---*/

//  Read_SU2_Restart_Metadata(geometry[MESH_0], config, true, restart_filename);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Read all lines in the restart file ---*/

  long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
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

  if (iPoint_Global_Local < nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars; Restart_Vars = nullptr;
  delete [] Restart_Data; Restart_Data = nullptr;

}
