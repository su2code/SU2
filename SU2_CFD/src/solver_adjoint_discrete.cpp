/*!
 * \file solver_adjoint_discrete.cpp
 * \brief Main subroutines for solving the discrete adjoint problem.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

CDiscAdjSolver::CDiscAdjSolver(void) : CSolver () {

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config)  : CSolver() {

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  unsigned short iVar, iMarker, iDim;

  bool restart = config->GetRestart();

  unsigned long iVertex, iPoint, index;
  string text_line, mesh_filename;
  ifstream restart_file;
  string filename, AdjExt;
  su2double dull_val;
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  /*--- Initialize arrays to NULL ---*/

  CSensitivity = NULL;

  Sens_Geo   = NULL;
  Sens_Mach  = NULL;
  Sens_AoA   = NULL;
  Sens_Press = NULL;
  Sens_Temp  = NULL;

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;


  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Allocate the node variables ---*/

  node = new CVariable*[nPoint];

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 1e-16;

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      CSensitivity[iMarker]        = new su2double [geometry->nVertex[iMarker]];
  }

  Sens_Geo  = new su2double[nMarker];
  Sens_Mach = new su2double[nMarker];
  Sens_AoA  = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Sens_Geo[iMarker]  = 0.0;
      Sens_Mach[iMarker] = 0.0;
      Sens_AoA[iMarker]  = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          CSensitivity[iMarker][iVertex] = 0.0;
      }
  }


  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  if (!restart || (iMesh != MESH_0)) {

    /*--- Restart the solution from zero ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);

  }
  else {

    /*--- Restart the solution from file information ---*/
    mesh_filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(mesh_filename);

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no file ---*/
    if (restart_file.fail()) {
      if (rank == MASTER_NODE) cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }

    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    
    map<unsigned long,unsigned long> Global2Local;
    map<unsigned long,unsigned long>::const_iterator MI;
    
    /*--- Now fill array with the transform values only for local points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }

    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- Skip coordinates ---*/
    unsigned short skipVars = nDim;

    /*--- Skip flow adjoint variables ---*/
    if (Kind_Solver == RUNTIME_TURB_SYS) {
      if (compressible) {
        skipVars += nDim + 2;
      }
      if (incompressible) {
        skipVars += nDim + 1;
      }
    }

    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
        
        iPoint_Local = Global2Local[iPoint_Global];
        
        point_line >> index;
        for (iVar = 0; iVar < skipVars; iVar++) { point_line >> dull_val;}
        for (iVar = 0; iVar < nVar; iVar++) { point_line >> Solution[iVar];}
        node[iPoint_Local] = new CDiscAdjVariable(Solution, nDim, nVar, config);
        iPoint_Global_Local++;
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

    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);
    }

    /*--- Close the restart file ---*/
    restart_file.close();

  }

  /*--- Store the direct solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetSolution_Direct(direct_solver->node[iPoint]->GetSolution());
  }
}

CDiscAdjSolver::~CDiscAdjSolver(void) { 

  unsigned short iMarker;

  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

  if (Sens_Geo   != NULL) delete [] Sens_Geo;
  if (Sens_Mach  != NULL) delete [] Sens_Mach;
  if (Sens_AoA   != NULL) delete [] Sens_AoA;
  if (Sens_Press != NULL) delete [] Sens_Press;
  if (Sens_Temp  != NULL) delete [] Sens_Temp;

}

void CDiscAdjSolver::SetRecording(CGeometry* geometry, CConfig *config, unsigned short kind_recording) {


  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->SetSolution(node[iPoint]->GetSolution_Direct());
  }

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
      }
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the meanflow iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjSolver::RegisterSolution(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, nPoint = geometry->GetnPoint();

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND,
  input = true;

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->RegisterSolution(input);
  }
  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      direct_solver->node[iPoint]->RegisterSolution_time_n();
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      direct_solver->node[iPoint]->RegisterSolution_time_n1();
    }
  }
}

void CDiscAdjSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) {

  /*--- Register farfield values as input ---*/

  if((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS)) {

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


    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/
}

void CDiscAdjSolver::RegisterOutput(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint, nPoint = geometry->GetnPoint();

  /*--- Register variables as output of the solver iteration ---*/

  bool input = false;

  /*--- Register output variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->RegisterSolution(input);
  }
}

void CDiscAdjSolver::RegisterObj_Func(CConfig *config) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Here we can add new (scalar) objective functions ---*/
  if (config->GetnObj()==1) {
    switch (config->GetKind_ObjFunc()) {
    case DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CD();
      if (config->GetFixed_CL_Mode()) ObjFunc_Value -= config->GetdCD_dCL() * direct_solver->GetTotal_CL();
      if (config->GetFixed_CM_Mode()) ObjFunc_Value -= config->GetdCD_dCM() * direct_solver->GetTotal_CMy();
      break;
    case LIFT_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CL();
      break;
    case AERO_DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_AeroCD();
      break;
    case RADIAL_DISTORTION:
      ObjFunc_Value = direct_solver->GetTotal_RadialDistortion();
      break;
    case CIRCUMFERENTIAL_DISTORTION:
      ObjFunc_Value = direct_solver->GetTotal_CircumferentialDistortion();
      break;
    case ELLIPTIC_SPANLOAD:
      ObjFunc_Value = direct_solver->GetTotal_EllipticDiff();
      break;
    case MAX_SECTIONAL_CL:
      ObjFunc_Value = direct_solver->GetTotal_MaxSecCL();
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
    case AVG_TOTAL_PRESSURE:
      ObjFunc_Value = direct_solver->GetOneD_TotalPress();
      break;
    case AVG_OUTLET_PRESSURE:
      ObjFunc_Value = direct_solver->GetOneD_FluxAvgPress();
      break;
    case MASS_FLOW_RATE:
      ObjFunc_Value = direct_solver->GetOneD_MassFlowRate();
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

void CDiscAdjSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config) {
  int rank = MASTER_NODE;

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

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) {
    SU2_TYPE::SetDerivative(ObjFunc_Value, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
  }
}

void CDiscAdjSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config) {

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  bool time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
      SetRes_RMS(iVar,0.0);
      SetRes_Max(iVar,0.0,0);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Set the old solution ---*/

    node[iPoint]->Set_OldSolution();

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution(Solution);
  }

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Extract the adjoint solution at time n ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n(Solution);

      /*--- Store the adjoint solution at time n ---*/

      node[iPoint]->Set_Solution_time_n(Solution);
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Extract the adjoint solution at time n-1 ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n1(Solution);

      /*--- Store the adjoint solution at time n-1 ---*/

      node[iPoint]->Set_Solution_time_n1(Solution);
    }
  }

  /*--- Set the residuals ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->GetSolution_Old(iVar);

          AddRes_RMS(iVar,residual*residual);
          AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_RMS(geometry, config);
}

void CDiscAdjSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) {
  su2double Local_Sens_Press, Local_Sens_Temp, Local_Sens_AoA, Local_Sens_Mach;

  /*--- Extract the adjoint values of the farfield values ---*/

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS)) {
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

  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/

}

void CDiscAdjSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config) {

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
      config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution[iVar] = node[iPoint]->GetSolution(iVar);
    }
    if (dual_time) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] += node[iPoint]->GetDual_Time_Derivative(iVar);
      }
    }
    direct_solver->node[iPoint]->SetAdjointSolution(Solution);
  }
}

void CDiscAdjSolver::SetSensitivity(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;
  su2double *Coord, Sensitivity, eps;

  bool time_stepping = (config->GetUnsteady_Simulation() != STEADY);

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Coord = geometry->node[iPoint]->GetCoord();

    for (iDim = 0; iDim < nDim; iDim++) {

      Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);

      /*--- Set the index manually to zero. ---*/

     AD::ResetInput(Coord[iDim]);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetLimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetSharpEdgesCoeff()*eps )
          Sensitivity = 0.0;
      }
      if (!time_stepping) {
        node[iPoint]->SetSensitivity(iDim, Sensitivity);
      } else {
        node[iPoint]->SetSensitivity(iDim, node[iPoint]->GetSensitivity(iDim) + Sensitivity);
      }
    }
  }
  SetSurface_Sensitivity(geometry, config);
}

void CDiscAdjSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker,iDim;
  unsigned long iVertex, iPoint;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area;
  su2double Total_Sens_Geo_local = 0.0;
  Total_Sens_Geo = 0.0;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker] = 0.0;
    /*--- Loop over boundary markers to select those for Euler walls and NS walls ---*/

    if(config->GetMarker_All_KindBC(iMarker) == EULER_WALL
       || config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
       || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Prod = 0.0;
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          /*--- retrieve the gradient calculated with AD -- */
          SensDim = node[iPoint]->GetSensitivity(iDim);

          /*--- calculate scalar product for projection onto the normal vector ---*/
          Prod += Normal[iDim]*SensDim;

          Area += Normal[iDim]*Normal[iDim];
        }

        Area = sqrt(Area);

        /*--- projection of the gradient
         *     calculated with AD onto the normal
         *     vector of the surface ---*/
        Sens = Prod/Area;

        /*--- Compute sensitivity for each surface point ---*/
        CSensitivity[iMarker][iVertex] = -Sens;
        if (geometry->node[iPoint]->GetDomain()) {
          Sens_Geo[iMarker] += Sens*Sens;
        }
      }
      Total_Sens_Geo_local += sqrt(Sens_Geo[iMarker]);
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Total_Sens_Geo_local,&Total_Sens_Geo,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#else
  Total_Sens_Geo = Total_Sens_Geo_local;
#endif
}

void CDiscAdjSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  bool dual_time_1st = (config_container->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  su2double *solution_n, *solution_n1;
  unsigned long iPoint;
  unsigned short iVar;
  if (dual_time) {
      for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++) {
          solution_n = node[iPoint]->GetSolution_time_n();
          solution_n1 = node[iPoint]->GetSolution_time_n1();

          for (iVar=0; iVar < nVar; iVar++) {
              node[iPoint]->SetDual_Time_Derivative(iVar, solution_n[iVar]+node[iPoint]->GetDual_Time_Derivative_n(iVar));
              node[iPoint]->SetDual_Time_Derivative_n(iVar, solution_n1[iVar]);

            }

        }

    }
}
