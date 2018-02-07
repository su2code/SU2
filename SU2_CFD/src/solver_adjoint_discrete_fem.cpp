/*!
 * \file solver_adjoint_discrete_fem.cpp
 * \brief Main subroutines for solving the finite element discrete adjoint problem.
 * \author B. Munguia
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CFEM_DG_DiscAdjSolver::CFEM_DG_DiscAdjSolver(void) : CSolver () {

}

CFEM_DG_DiscAdjSolver::CFEM_DG_DiscAdjSolver(CGeometry *geometry, CConfig *config)  : CSolver() {

}

CFEM_DG_DiscAdjSolver::CFEM_DG_DiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  unsigned short iVar, iMarker, iDim;
  unsigned long iDOF;
  string text_line, mesh_filename;
  ifstream restart_file;
  string filename, AdjExt;

  bool fsi = config->GetFSI_Simulation();

  nVar    = direct_solver->GetnVar();
  nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  /*--- Initialize arrays to NULL ---*/

  CSensitivity = NULL;

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  nDOFsLocOwned = geometry->GetnPointDomain();
  nDOFsGlobal   = geometry->GetGlobal_nPointDomain();

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  nMeshPoints = DGGeometry->GetNMeshPoints();
  meshPoints  = DGGeometry->GetMeshPoints();

  nDOFsLocTot = nDOFsLocOwned;
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) nDOFsLocTot += volElem[i].nDOFsSol;

  boundaries = DGGeometry->GetBoundaries();

  nStandardBoundaryFacesSol  = DGGeometry->GetNStandardBoundaryFacesSol();
  standardBoundaryFacesSol  = DGGeometry->GetStandardBoundaryFacesSol();

  nDOFsBoundary = new unsigned long[nMarker];

  for(iMarker = 0; iMarker < nMarker; iMarker++){
    const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
    const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

    nDOFsBoundary[iMarker] = 0;

    /*--- Loop over the faces of this boundary. ---*/
    for(unsigned long l=0; l<nSurfElem; ++l) {

      /* Get the required information from the corresponding standard face. */
      const unsigned short ind       = surfElem[l].indStandardElement;
      const unsigned short nDOFsFace = standardBoundaryFacesSol[ind].GetNDOFsFace();
      nDOFsBoundary[iMarker] += nDOFsFace;
    }
  }

  /*--- Allocate solution vectors ---*/

  VecSolDOFs.resize(nVar*nDOFsLocOwned);
  VecSolDOFsNew.resize(nVar*nDOFsLocOwned);
  VecSolDOFsDirect.resize(nVar*nDOFsLocOwned);
  VecSolDOFsSens.resize(nDim*nDOFsLocOwned);

  unsigned long ii = 0;
  for (iDOF = 0; iDOF < nDOFsLocOwned; iDOF++){
    ii = nVar*iDOF;
    for (iVar = 0; iVar < nVar; iVar++, ii++){
      VecSolDOFs[ii+iVar] = 1e-16;
      VecSolDOFsNew[ii+iVar] = 1e-16;
    }
    for(iDim = 0; iDim < nDim; iDim++){
      VecSolDOFsSens[nDim*iDOF+iDim] = 1e-16;
    }
  }

  StoreSolution_Direct();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar];            for(iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 1e-16;

  /*--- Define some auxiliary vectors related to the geometry adjoint (nDim) ---*/
  Solution_Geometry = new su2double[nDim];     for (iDim = 0; iDim < nDim; iDim++) Solution_Geometry[iDim] = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      CSensitivity[iMarker]        = new su2double [nDOFsBoundary[iMarker]];
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iDOF = 0; iDOF < nDOFsBoundary[iMarker]; iDOF++) {
          CSensitivity[iMarker][iDOF] = 0.0;
      }
  }

}

CFEM_DG_DiscAdjSolver::~CFEM_DG_DiscAdjSolver(void) { 

  unsigned short iMarker;

  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

}

void CFEM_DG_DiscAdjSolver::StoreSolution_Direct(){

  su2double *Solution = direct_solver->GetVecSolDOFs();

  for(unsigned long iDOF = 0; iDOF < nDOFsLocOwned; iDOF++){
    unsigned long ii = nVar*iDOF;
    for(unsigned short iVar = 0; iVar < nVar; iVar++){
      VecSolDOFsDirect[ii+iVar] = Solution[ii+iVar];
    }
  }
}

void CFEM_DG_DiscAdjSolver::SetRecording(CGeometry* geometry, CConfig *config){


  bool time_stepping  = config->GetUnsteady_Simulation() == TIME_STEPPING;

  unsigned long iDOF;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  direct_solver->ResetSolution_Direct(VecSolDOFsDirect);

  if (time_stepping) {
    for (iDOF = 0; iDOF < nDOFsLocOwned; iDOF++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        //TODO: reset input for unsteady
        //AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CFEM_DG_DiscAdjSolver::SetMesh_Recording(CGeometry** geometry, CVolumetricMovement *grid_movement, CConfig *config) {


//  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
//      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
//  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

//  unsigned long ExtIter = config->GetExtIter();

  unsigned long iPoint;
  unsigned short iDim;

  /*--- Reset the solution to the initial (converged) position ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iDim = 0; iDim < nDim; iDim++){
      geometry[MESH_0]->node[iPoint]->SetCoord(iDim,node[iPoint]->GetGeometry_Direct(iDim));
    }
  }

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  grid_movement->UpdateDualGrid(geometry[MESH_0], config);

  /*--- After updating the dual mesh, compute the grid velocities (only dynamic problems). ---*/
//  if (time_n_needed){
//    geometry[MESH_0]->SetGridVelocity(config, ExtIter);
//  }

  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels. ---*/

  grid_movement->UpdateMultiGrid(geometry, config);

//  if (time_n_needed){
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//      for (iVar = 0; iVar < nVar; iVar++){
//        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
//      }
//    }
//  }
//  if (time_n1_needed){
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//      for (iVar = 0; iVar < nVar; iVar++){
//        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
//      }
//    }
//  }

}

void CFEM_DG_DiscAdjSolver::RegisterSolution(CGeometry *geometry, CConfig *config) {
  
  unsigned long iDOF, ii;
  unsigned short iVar;

  bool time_stepping  = config->GetUnsteady_Simulation() == TIME_STEPPING;

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  for (iDOF = 0; iDOF < nDOFsLocOwned; iDOF++) {
    direct_solver->RegisterSolution(iDOF,true);
  }
  if (time_stepping) {
    for (iDOF = 0; iDOF < nDOFsLocOwned; iDOF++) {
      //direct_solver->node[iPoint]->RegisterSolution_time_n();
    }
  }
}

void CFEM_DG_DiscAdjSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) {

  /*--- Register farfield values as input ---*/

  if((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS && !config->GetBoolTurbomachinery())) {

    su2double Velocity_Ref = config->GetVelocity_Ref();
    Alpha                  = config->GetAoA()*PI_NUMBER/180.0;
    Beta                   = config->GetAoS()*PI_NUMBER/180.0;
    Mach                   = config->GetMach();
    Pressure               = config->GetPressure_FreeStreamND();
    Temperature            = config->GetTemperature_FreeStreamND();

    su2double SoundSpeed = 0.0;
    
    if (nDim == 2) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*Mach); }
    if (nDim == 3) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*cos(Beta)*Mach); }

    if (!reset) {
      AD::RegisterInput(Mach);
      AD::RegisterInput(Alpha);
      AD::RegisterInput(Temperature);
      AD::RegisterInput(Pressure);
    }

    /*--- Recompute the free stream velocity ---*/

    if (nDim == 2) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Alpha)*Mach*SoundSpeed/Velocity_Ref;
    }
    if (nDim == 3) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[2] = sin(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
    }

    config->SetTemperature_FreeStreamND(Temperature);
    direct_solver->SetTemperature_Inf(Temperature);
    config->SetPressure_FreeStreamND(Pressure);
    direct_solver->SetPressure_Inf(Pressure);

  }

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && config->GetBoolTurbomachinery()){

    BPressure = config->GetPressureOut_BC();
    Temperature = config->GetTotalTemperatureIn_BC();

    if (!reset){
      AD::RegisterInput(BPressure);
      AD::RegisterInput(Temperature);
    }

    config->SetPressureOut_BC(BPressure);
    config->SetTotalTemperatureIn_BC(Temperature);
  }


    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/
}

void CFEM_DG_DiscAdjSolver::RegisterOutput(CGeometry *geometry, CConfig *config) {

  unsigned long iDOF, ii;
  unsigned short iVar;

  bool time_stepping  = config->GetUnsteady_Simulation() == TIME_STEPPING;

  /*--- Register output variables on the tape ---*/

  for (iDOF = 0; iDOF < nDOFsLocOwned; iDOF++) {
    direct_solver->RegisterSolution(iDOF,false);
  }
  if (time_stepping) {
    for (iDOF = 0; iDOF < nDOFsLocOwned; iDOF++) {
      //direct_solver->node[iPoint]->RegisterSolution_time_n();
    }
  }
}

void CFEM_DG_DiscAdjSolver::RegisterObj_Func(CConfig *config) {

  /*--- Here we can add new (scalar) objective functions ---*/
  if (config->GetnObj()==1) {
    switch (config->GetKind_ObjFunc()) {
    case DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CD();
      if (config->GetFixed_CL_Mode()) ObjFunc_Value -= config->GetdCD_dCL() * direct_solver->GetTotal_CL();
      if (config->GetFixed_CM_Mode()) ObjFunc_Value -= config->GetdCD_dCMy() * direct_solver->GetTotal_CMy();
      break;
    case LIFT_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CL();
      break;
    case SIDEFORCE_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CSF();
      break;
    case EFFICIENCY:
      ObjFunc_Value = direct_solver->GetTotal_CEff();
      break;
    case MOMENT_X_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMx();
      break;
    case MOMENT_Y_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMy();
      break;
    case MOMENT_Z_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMz();
      break;
    case EQUIVALENT_AREA:
      ObjFunc_Value = direct_solver->GetTotal_CEquivArea();
      break;
    }

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
  else{
    ObjFunc_Value = direct_solver->GetTotal_ComboObj();
  }
  if (rank == MASTER_NODE) {
    AD::RegisterOutput(ObjFunc_Value);
  }
}

void CFEM_DG_DiscAdjSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config) {

  bool time_stepping = config->GetUnsteady_Simulation() != STEADY;
  unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
  unsigned long ExtIter = config->GetExtIter();
  su2double seeding = 1.0;

  if (time_stepping) {
    if (ExtIter < IterAvg_Obj) {
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else {
      seeding = 0.0;
    }
  }

  if (rank == MASTER_NODE) {
    SU2_TYPE::SetDerivative(ObjFunc_Value, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
  }
}

void CFEM_DG_DiscAdjSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  bool time_stepping  = config->GetUnsteady_Simulation() == TIME_STEPPING;

  unsigned short iVar;
  unsigned long iDOF, ii;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
      SetRes_RMS(iVar,0.0);
      SetRes_Max(iVar,0.0,0);
  }


  /*--- Set the old solution ---*/

  memcpy(VecSolDOFs.data(), VecSolDOFsNew.data(), VecSolDOFs.size()*sizeof(su2double));

  /*--- Extract the adjoint solution ---*/

  direct_solver->GetAdjointSolution(VecSolDOFsNew);

  if (time_stepping) {
    // TODO: code for unsteady
  }

  /*--- Set the residuals ---*/

  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Store the coordinate of the first vertex of this element to give an
       indication for the location of the maximum residual. */
    const unsigned long ind = volElem[l].nodeIDsGrid[0];
    const su2double *coor   = meshPoints[ind].coor;

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *solDOFsOld = VecSolDOFs.data() + offset;
    const su2double *solDOFsNew = VecSolDOFsNew.data() + offset;

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      for(unsigned short iVar=0; iVar<nVar; ++iVar, ++i) {
        residual = solDOFsNew[i] - solDOFsOld[i];

        AddRes_RMS(iVar, residual*residual);
        AddRes_Max(iVar, fabs(residual), globalIndex, coor);
      }
    }
  }

  /*--- Compute the root mean square residual. Note that the SetResidual_RMS
        function cannot be used, because that is for the FV solver.    ---*/

#ifdef HAVE_MPI
  /*--- Parallel mode. The local L2 norms must be added to obtain the
        global value. Also check for divergence. ---*/
  vector<su2double> rbuf(nVar);

  /*--- Disable the reduce for the residual to avoid overhead if requested. ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(Residual_RMS, rbuf.data(), nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {

      if (rbuf[iVar] != rbuf[iVar])
        SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

      SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf[iVar]/nDOFsGlobal)));
    }
  }

#else
  /*--- Sequential mode. Check for a divergence of the solver and compute
        the L2-norm of the residuals. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if(GetRes_RMS(iVar) != GetRes_RMS(iVar))
      SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(GetRes_RMS(iVar)/nDOFsGlobal)));
  }

#endif

}

void CFEM_DG_DiscAdjSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) {

  /*--- Extract the adjoint values of the farfield values ---*/

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && !config->GetBoolTurbomachinery()) {
    su2double Local_Sens_Press, Local_Sens_Temp, Local_Sens_AoA, Local_Sens_Mach;

    Local_Sens_Mach  = SU2_TYPE::GetDerivative(Mach);
    Local_Sens_AoA   = SU2_TYPE::GetDerivative(Alpha);
    Local_Sens_Temp  = SU2_TYPE::GetDerivative(Temperature);
    Local_Sens_Press = SU2_TYPE::GetDerivative(Pressure);

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&Local_Sens_Mach,  &Total_Sens_Mach,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_AoA,   &Total_Sens_AoA,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Temp,  &Total_Sens_Temp,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    Total_Sens_Mach  = Local_Sens_Mach;
    Total_Sens_AoA   = Local_Sens_AoA;
    Total_Sens_Temp  = Local_Sens_Temp;
    Total_Sens_Press = Local_Sens_Press;
#endif
  }

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && config->GetBoolTurbomachinery()){
    su2double Local_Sens_BPress, Local_Sens_Temperature;

    Local_Sens_BPress = SU2_TYPE::GetDerivative(BPressure);
    Local_Sens_Temperature = SU2_TYPE::GetDerivative(Temperature);

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&Local_Sens_BPress,   &Total_Sens_BPress,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Temperature,   &Total_Sens_Temp,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    Total_Sens_BPress = Local_Sens_BPress;
    Total_Sens_Temp = Local_Sens_Temperature;
#endif

  }

  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/

}

void CFEM_DG_DiscAdjSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config) {


  direct_solver->SetAdjointSolution(VecSolDOFsNew);

}

void CFEM_DG_DiscAdjSolver::SetSensitivity(CGeometry *geometry, CConfig *config) {

  unsigned long iDOF;
  unsigned short iDim;
  su2double *Coord, Sensitivity, eps;

  bool time_stepping = (config->GetUnsteady_Simulation() != STEADY);

  /* Loop over the owned elements. */
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Loop over the DOFs of this element. */
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

      // Set the pointer to the solution of this DOF and to the
      // coordinates of its corresponding node ID of the grid.
      const unsigned long ind = volElem[i].nodeIDsGrid[j];
      su2double *coor   = meshPoints[ind].coor;

      for(iDim = 0; iDim < nDim; iDim++){
        if(!time_stepping){
          VecSolDOFsSens[ind*nDim+iDim] = SU2_TYPE::GetDerivative(coor[iDim]);
        }
        else{
          VecSolDOFsSens[ind*nDim+iDim] += SU2_TYPE::GetDerivative(coor[iDim]);
        }

        /*--- Set the index manually to zero. ---*/

        AD::ResetInput(coor[iDim]);
      }
    }
  }
  SetSurface_Sensitivity(geometry, config);
}

void CFEM_DG_DiscAdjSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, iDim, iMarker_Monitoring;
  unsigned long iVertex, iPoint;
  su2double Prod, Sens = 0.0, SensDim, Area, Sens_Vertex, *Sens_Geo;
  Total_Sens_Geo = 0.0;
  string Monitoring_Tag, Marker_Tag;

  Sens_Geo = new su2double[config->GetnMarker_Monitoring()];
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Sens_Geo[iMarker_Monitoring] = 0.0;
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- Loop over boundary markers to select those for Euler walls and NS walls ---*/

    if(config->GetMarker_All_KindBC(iMarker) == EULER_WALL
       || config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
       || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {

      Sens = 0.0;
      iVertex = 0;

      /* Easier storage of the boundary faces for this boundary marker. */
      const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
      const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

      /*--- Loop over the faces of this boundary. ---*/
      for(unsigned long l=0; l<nSurfElem; ++l) {

        /* Get the required information from the corresponding standard face. */
        const unsigned short ind   = surfElem[l].indStandardElement;
        const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
        const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();

        /* Loop over the integration points of this surface element. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the solution, the normals and the coordinates
             for this integration point. */
          const su2double *normals = surfElem[l].metricNormalsFace.data()
                                   + i*(nDim+1);
          const su2double *Coord   = surfElem[l].coorIntegrationPoints.data()
                                   + i*nDim;  

          Prod = 0.0;
          Area = 0.0;
          const unsigned long  volID     = surfElem[l].volElemID;
          for (iDim = 0; iDim < nDim; iDim++) {
            /*--- retrieve the gradient calculated with AD -- */
            SensDim = VecSolDOFsSens[volID*nDim+iDim];

            /*--- calculate scalar product for projection onto the normal vector ---*/
            Prod += normals[iDim]*SensDim;

            Area += normals[iDim]*normals[iDim];
          }

          Area = sqrt(Area);

          /*--- Projection of the gradient calculated with AD onto the normal vector of the surface ---*/

          Sens_Vertex = Prod/Area;
          CSensitivity[iMarker][iVertex] = -Sens_Vertex;
          Sens += Sens_Vertex*Sens_Vertex;

          iVertex++;
        }
      }

      if (config->GetMarker_All_Monitoring(iMarker) == YES){

        /*--- Compute sensitivity for each surface point ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Sens_Geo[iMarker_Monitoring] = Sens;
          }
        }
      }
    }
  }

#ifdef HAVE_MPI
  su2double *MySens_Geo;
  MySens_Geo = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySens_Geo[iMarker_Monitoring] = Sens_Geo[iMarker_Monitoring];
    Sens_Geo[iMarker_Monitoring]   = 0.0;
  }

  SU2_MPI::Allreduce(MySens_Geo, Sens_Geo, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  delete [] MySens_Geo;
#endif

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Sens_Geo[iMarker_Monitoring] = sqrt(Sens_Geo[iMarker_Monitoring]);
    Total_Sens_Geo   += Sens_Geo[iMarker_Monitoring];
  }
  
  delete [] Sens_Geo;

}

void CFEM_DG_DiscAdjSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  bool time_stepping = (config_container->GetUnsteady_Simulation() == TIME_STEPPING);
  su2double *solution_n, *solution_n1;
  unsigned long iDOF;
  unsigned short iVar;
  if (time_stepping) {
    // TODO: preprocess for unsteady  
  }
}

void CFEM_DG_DiscAdjSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iVar;
  unsigned long index;

  string UnstExt, text_line;
  ifstream restart_file;

  string restart_filename, filename;

  const bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  const bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  filename = config->GetSolution_AdjFileName();
  restart_filename = config->GetObjFunc_Extension(filename);

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned short rbuf_NotMatching = 0;
  unsigned long nDOF_Read = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

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
      for (iVar = 0; iVar < nVar; iVar++) {
        VecSolDOFs[nVar*iPoint_Local+iVar] = Restart_Data[index+iVar];
      }
      /*--- Update the local counter nDOF_Read. ---*/
      ++nDOF_Read;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/
  if(nDOF_Read < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
  unsigned short sbuf_NotMatching = rbuf_NotMatching;
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, MPI_COMM_WORLD);
#endif

  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename.data() +
                   string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."),
                   CURRENT_FUNCTION);

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}
