/*!
 * \file CDiscAdjFEASolver.cpp
 * \brief Main subroutines for solving adjoint FEM elasticity problems.
 * \author R. Sanchez
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

CDiscAdjFEASolver::CDiscAdjFEASolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver,
                                     unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  adjoint = true;

  unsigned long iPoint;

  const bool dynamic = (config->GetTime_Domain());

  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS.resize(nVar,1.0);
  Residual_Max.resize(nVar,1.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()) {

    Residual_BGS.resize(nVar,1.0);
    Residual_Max_BGS.resize(nVar,1.0);
    Point_Max_BGS.resize(nVar,0);
    Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
  }

  /*--- Initialize the adjoint solution. ---*/

  vector<su2double> init(nVar,1e-16);
  nodes = new CDiscAdjFEABoundVariable(init.data(), init.data(), init.data(), nPoint, nDim, nVar, dynamic, config);
  SetBaseClassPointerToNodes();

  /*--- Set which points are vertices and allocate boundary data. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      long iVertex = geometry->nodes->GetVertex(iPoint, iMarker);
      if (iVertex >= 0) {
        nodes->Set_isVertex(iPoint,true);
        break;
      }
    }
  nodes->AllocateBoundaryVariables(config);

  /*--- Initialize vector structures for multiple material definition ---*/

  nMPROP = config->GetnElasticityMod();

  /*--- For a material to be fully defined, we need to have the same number for all three parameters ---*/
  bool checkDef = ((config->GetnElasticityMod() == config->GetnPoissonRatio()) &&
                   (config->GetnElasticityMod() == config->GetnMaterialDensity()) &&
                   (config->GetnMaterialDensity() == config->GetnPoissonRatio()));

  if (!checkDef){
    SU2_MPI::Error("WARNING: For a material to be fully defined, E, Nu and Rho need to have the same dimensions.", CURRENT_FUNCTION);
  }

  E.resize(nMPROP);
  Nu.resize(nMPROP);
  Rho.resize(nMPROP);
  Rho_DL.resize(nMPROP);

  /*--- Initialize vector structures for multiple electric regions ---*/

  de_effects = config->GetDE_Effects();

  if (de_effects) {
    nEField = config->GetnElectric_Field();
    EField.resize(nEField);
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

  if (fea_dv) ReadDV(config);

}

CDiscAdjFEASolver::~CDiscAdjFEASolver() { delete nodes; }

void CDiscAdjFEASolver::SetRecording(CGeometry* geometry, CConfig *config){

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint, nodes->GetSolution_Direct(iPoint));
  }

  if (config->GetTime_Domain()){
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

  const bool input = true;
  const bool dynamic = config->GetTime_Domain();
  const bool push_index = !config->GetMultizone_Problem();

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

    const bool pseudo_static = config->GetPseudoStatic();

    for (iVar = 0; iVar < nMPROP; iVar++) {
      E[iVar]      = config->GetElasticyMod(iVar);
      Nu[iVar]     = config->GetPoissonRatio(iVar);
      Rho[iVar]    = pseudo_static? 0.0 : config->GetMaterialDensity(iVar);
      Rho_DL[iVar] = config->GetMaterialDensity(iVar);
    }

    /*--- Read the values of the electric field ---*/
    if (de_effects) {
      for (iVar = 0; iVar < nEField; iVar++)
        EField[iVar] = config->Get_Electric_Field_Mod(iVar);
    }

    /*--- Reset index, otherwise messes up other derivatives ---*/
    if (fea_dv) {
      for (iVar = 0; iVar < nDV; iVar++) AD::ResetInput(DV[iVar]);
    }

    if (!reset) {
      const bool local_index = config->GetMultizone_Problem();
      const bool push_index = !local_index;

      E.Register(push_index);
      Nu.Register(push_index);
      Rho.Register(push_index);
      Rho_DL.Register(push_index);
      if (de_effects) EField.Register(push_index);
      if (fea_dv) DV.Register(push_index);

      /*--- Explicitly store the tape indices for when we extract the derivatives ---*/
      if (local_index) {
        E.SetIndex();
        Nu.SetIndex();
        Rho.SetIndex();
        Rho_DL.SetIndex();
        if (de_effects) EField.SetIndex();
        if (fea_dv) DV.SetIndex();
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

  const bool input = false;
  const bool dynamic = config->GetTime_Domain();
  const bool push_index = !config->GetMultizone_Problem();

  /*--- Register variables as output of the solver iteration ---*/

  direct_solver->GetNodes()->RegisterSolution(input, push_index);

  if (dynamic) {
    /*--- Register acceleration (u'') and velocity (u') at time step n ---*/
    direct_solver->GetNodes()->RegisterSolution_Accel(input);
    direct_solver->GetNodes()->RegisterSolution_Vel(input);
  }

}

void CDiscAdjFEASolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  const bool dynamic = config->GetTime_Domain();
  const bool multizone = config->GetMultizone_Problem();

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  su2double Solution[MAXNVAR] = {0.0}, Solution_Vel[MAXNVAR] = {0.0}, Solution_Accel[MAXNVAR] = {0.0};

  /*--- Set Residuals to zero ---*/

  SetResToZero();

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

      Residual_RMS[iVar] += residual*residual;
      AddRes_Max(iVar,fabs(residual),geometry->nodes->GetGlobalIndex(iPoint),geometry->nodes->GetCoord(iPoint));
    }
    if (dynamic){
      for (iVar = 0; iVar < nVar; iVar++){
        residual = nodes->GetSolution_Accel(iPoint, iVar) - nodes->GetSolution_Old_Accel(iPoint, iVar);

        Residual_RMS[iVar] += residual*residual;
        AddRes_Max(iVar,fabs(residual),geometry->nodes->GetGlobalIndex(iPoint),geometry->nodes->GetCoord(iPoint));
      }
      for (iVar = 0; iVar < nVar; iVar++){
        residual = nodes->GetSolution_Vel(iPoint, iVar) - nodes->GetSolution_Old_Vel(iPoint, iVar);

        Residual_RMS[iVar] += residual*residual;
        AddRes_Max(iVar,fabs(residual),geometry->nodes->GetGlobalIndex(iPoint),geometry->nodes->GetCoord(iPoint));
      }
    }
  }

  SetResidual_RMS(geometry, config);

  SetIterLinSolver(direct_solver->System.GetIterations());
  SetResLinSolver(direct_solver->System.GetResidual());

}

void CDiscAdjFEASolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){

  if (KindDirect_Solver != RUNTIME_FEA_SYS) return;

  /*--- Sensitivities of material properties and design variables. ---*/

  E.GetDerivative();
  Nu.GetDerivative();
  Rho.GetDerivative();
  Rho_DL.GetDerivative();
  if (de_effects) EField.GetDerivative();
  if (fea_dv) DV.GetDerivative();

  /*--- Extract the flow traction sensitivities. ---*/

  if (config->GetnMarker_Fluid_Load() > 0) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
      for (unsigned short iDim = 0; iDim < nDim; iDim++){
        su2double val_sens = direct_solver->GetNodes()->ExtractFlowTraction_Sensitivity(iPoint,iDim);
        nodes->SetFlowTractionSensitivity(iPoint, iDim, val_sens);
      }
    }
  }

}

void CDiscAdjFEASolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config){

  const bool dynamic = (config->GetTime_Domain());
  const bool deform_mesh = (config->GetnMarker_Deform_Mesh() > 0);

  su2double Solution[MAXNVAR] = {0.0}, Solution_Vel[MAXNVAR] = {0.0}, Solution_Accel[MAXNVAR] = {0.0};

  unsigned short iVar;

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++){
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
        Solution[iVar] += nodes->GetDynamic_Derivative_n(iPoint,iVar);
      }
      for (iVar = 0; iVar < nVar; iVar++){
        Solution_Accel[iVar] = nodes->GetSolution_Accel(iPoint,iVar) + nodes->GetDynamic_Derivative_Accel_n(iPoint,iVar);
      }
      for (iVar = 0; iVar < nVar; iVar++){
        Solution_Vel[iVar] = nodes->GetSolution_Vel(iPoint,iVar) + nodes->GetDynamic_Derivative_Vel_n(iPoint,iVar);
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

  unsigned short iVar;

  if (config_container->GetTime_Domain()){
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++){
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

void CDiscAdjFEASolver::SetSensitivity(CGeometry *geometry, CConfig *config, CSolver*){

  E.UpdateTotal();
  Nu.UpdateTotal();
  Rho.UpdateTotal();
  Rho_DL.UpdateTotal();
  if (de_effects) EField.UpdateTotal();
  if (fea_dv) DV.UpdateTotal();

  /*--- Extract the topology optimization density sensitivities. ---*/

  direct_solver->ExtractAdjoint_Variables(geometry, config);

  /*--- Extract the geometric sensitivities ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    su2double *Coord = geometry->nodes->GetCoord(iPoint);

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {

      su2double Sensitivity;

      if(config->GetMultizone_Problem()) {
        Sensitivity = geometry->nodes->GetAdjointSolution(iPoint, iDim);
      }
      else {
        Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);
        /*--- Set the index manually to zero. ---*/
        AD::ResetInput(Coord[iDim]);
      }

      nodes->SetSensitivity(iPoint, iDim, Sensitivity);
    }
  }
}

void CDiscAdjFEASolver::ReadDV(const CConfig *config) {

  string filename;
  ifstream properties_file;

  /*--- Choose the filename of the design variable ---*/

  switch (config->GetDV_FEA()) {
    case YOUNG_MODULUS:
      filename = "dv_young.opt";
      break;
    case POISSON_RATIO:
      filename = "dv_poisson.opt";
      break;
    case DENSITY_VAL:
    case DEAD_WEIGHT:
      filename = "dv_density.opt";
      break;
    case ELECTRIC_FIELD:
      filename = "dv_efield.opt";
      break;
    default:
      filename = "dv.opt";
      break;
  }

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {

    if (rank == MASTER_NODE)
      cout << "There is no design variable file." << endl;

    nDV = 1;
    DV.resize(nDV);
    DV[0] = 1.0;
  }
  else{

    string text_line;
    vector<su2double> values;

    /*--- Skip the first line: it is the header ---*/
    getline (properties_file, text_line);

    while (getline (properties_file, text_line)) {
      istringstream point_line(text_line);

      unsigned long index;
      su2double value;
      point_line >> index >> value;

      values.push_back(value);
    }

    nDV = values.size();
    DV.resize(nDV);
    unsigned short iDV = 0;
    for (auto x : values) DV[iDV++] = x;

  }

}

void CDiscAdjFEASolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  auto filename = config->GetSolution_AdjFileName();
  auto restart_filename = config->GetObjFunc_Extension(filename);
  restart_filename = config->GetFilename(restart_filename, "", val_iter);

  BasicLoadRestart(geometry[MESH_0], config, restart_filename, geometry[MESH_0]->GetnDim());

}
