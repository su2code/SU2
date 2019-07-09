/*!
 * \file solution_direct_transition.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author A. Aranake
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
#include "../include/variables/CTransLMVariable.hpp"
#include "../include/variables/CTurbSAVariable.hpp"

CTransLMSolver::CTransLMSolver(void) : CTurbSolver() {}

CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint, index;
  su2double Density_Inf, Viscosity_Inf, tu_Inf, nu_tilde_Inf, Factor_nu_Inf, dull_val, rey;
  ifstream restart_file;
  char *cstr;
  string text_line;
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  
  cout << "Entered constructor for CTransLMSolver -AA\n";
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constans in the solver structure ---*/
  nDim = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  node = new CVariable*[geometry->GetnPoint()];
  
  /*--- Dimension of the problem --> 2 Transport equations (intermittency, Reth) ---*/
  nVar = 2;
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  if (iMesh == MESH_0) {
    
    /*--- Define some auxillary vectors related to the residual ---*/
    Residual     = new su2double[nVar]; Residual_RMS = new su2double[nVar];
    Residual_i   = new su2double[nVar]; Residual_j   = new su2double[nVar];
    Residual_Max = new su2double[nVar];

    /*--- Define some structures for locating max residuals ---*/
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }

    /*--- Define some auxiliar vector related with the solution ---*/
    Solution   = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
        
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {

      /*--- Point to point Jacobians ---*/
      Jacobian_i = new su2double* [nVar];
      Jacobian_j = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar] = new su2double [nVar];
        Jacobian_j[iVar] = new su2double [nVar];
      }
      /*--- Initialization of the structure of the whole Jacobian ---*/
      Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
      
      if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
          (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
        nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
        if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
      }
      
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
  
  /*--- Read farfield conditions from config ---*/
  Density_Inf       = config->GetDensity_FreeStreamND();
  Viscosity_Inf     = config->GetViscosity_FreeStreamND();
  Intermittency_Inf = config->GetIntermittency_FreeStream();
  tu_Inf            = config->GetTurbulenceIntensity_FreeStream();
  
  /*-- Initialize REth from correlation --*/
  if (tu_Inf <= 1.3) {
    REth_Inf = (1173.51-589.428*tu_Inf+0.2196/(tu_Inf*tu_Inf));
  } else {
    REth_Inf = 331.5*pow(tu_Inf-0.5658,-0.671);
  }
  rey = config->GetReynolds();

//  REth_Inf *= mach/rey;
  cout << "REth_Inf = " << REth_Inf << ", rey: "<< rey << " -AA" << endl;
  
  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/
  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
    
  /*--- Restart the solution from file information ---*/
  if (!restart) {
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      // TODO: Erase this bubble of specially initialized points -AA
//      if (iPoint == 9745||iPoint == 9746||iPoint == 9608||iPoint == 9609) {
//        node[iPoint] = new CTransLMVariable(nu_tilde_Inf, 0.0, 1100.0, nDim, nVar, config);
//      } else {
        node[iPoint] = new CTransLMVariable(nu_tilde_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);
 //     }
    }
    }

  }
  else {
    cout << "No LM restart yet!!" << endl; // TODO, Aniket
    int j;
    cin >> j;
    string mesh_filename = config->GetSolution_FlowFileName();
    cstr = new char [mesh_filename.size()+1];
    strcpy (cstr, mesh_filename.c_str());
    restart_file.open(cstr, ios::in);
    if (restart_file.fail()) {
      SU2_MPI::Error("There is no turbulent restart file.", CURRENT_FUNCTION);
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      getline(restart_file, text_line);
      istringstream point_line(text_line);
      if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
      if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
      node[iPoint] = new CTurbSAVariable(Solution[0], 0, nDim, nVar, config);
    }
    restart_file.close();
  }

}

CTransLMSolver::~CTransLMSolver(void) {
  
}

void CTransLMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
    LinSysRes.SetBlock_Zero(iPoint);
  Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CTransLMSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  /*--- Correction for separation-induced transition, Replace intermittency with gamma_eff ---*/
  for (unsigned int iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
    node[iPoint]->SetGammaEff();
}

void CTransLMSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Delta_flow, Vol;
  
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta_flow = Vol / (solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Delta = Delta_flow;
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      
      /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
      
      LinSysRes[total_index] = -LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]*Vol);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++)
      node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
  }
  
  /*--- MPI solution ---*/
    
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CTransLMSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  su2double *trans_var_i, *trans_var_j, *U_i, *U_j;
  unsigned long iEdge, iPoint, jPoint;

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge and normal vectors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
    U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
    numerics->SetConservative(U_i, U_j);

    /*--- Transition variables w/o reconstruction ---*/
    trans_var_i = node[iPoint]->GetSolution();
    trans_var_j = node[jPoint]->GetSolution();
    numerics->SetTransVar(trans_var_i, trans_var_j);

    /*--- Add and subtract Residual ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);

    /*--- Implicit part ---*/
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);

  }

}


void CTransLMSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());
    
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                            solver_container[FLOW_SOL]->node[jPoint]->GetSolution());
    
    /*--- Laminar Viscosity ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    /*--- Eddy Viscosity ---*/
    numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                             solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());
    
    /*--- Transition variables w/o reconstruction, and its gradients ---*/
    numerics->SetTransVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTransVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    numerics->SetConsVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient(),
                               solver_container[FLOW_SOL]->node[jPoint]->GetGradient());
    
    
    /*--- Compute residual, and Jacobians ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
}

void CTransLMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                       CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  su2double gamma_sep = 0.0;

  //cout << "Setting Trans residual -AA " << endl;
  cout << "\nBeginAA" << endl;
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    cout << "\niPoint: " << iPoint << endl;
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Laminar and eddy viscosity ---*/
    
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
    numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);
    
    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetTransVar(node[iPoint]->GetSolution(), NULL);
    // numerics->SetTransVarGradient(node[iPoint]->GetGradient(), NULL);  // Is this needed??
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set distance to the surface ---*/
    
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual_TransLM(Residual, Jacobian_i, NULL, config, gamma_sep);
    
    /*-- Store gamma_sep in variable class --*/
    
    node[iPoint]->SetGammaSep(gamma_sep);

    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }
}

void CTransLMSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh) {
}

void CTransLMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *U_i;
  su2double *U_domain = new su2double[nVar];
  su2double *U_wall   = new su2double[nVar];
  su2double *Normal   = new su2double[nDim];
  su2double *Residual = new su2double[nVar];

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT); 

//  cout << "Setting wall BC -AA\n";
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Set both interior and exterior point to current value ---*/
      for (iVar=0; iVar < nVar; iVar++) {
        U_domain[iVar] = node[iPoint]->GetSolution(iVar);
        U_wall[iVar]   = node[iPoint]->GetSolution(iVar);   
      }

      /*--- Set various quantities in the solver class ---*/
      numerics->SetNormal(Normal);
      numerics->SetTransVar(U_domain,U_wall);
      U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
      numerics->SetConservative(U_i, U_i);

      /*--- Compute the residual using an upwind scheme ---*/
//      cout << "BC calling SetResidual: -AA" << endl;
      numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//      cout << "Returned from BC call of SetResidual: -AA" << endl;
//      cout << "Residual[0] = " << Residual[0] << endl;
//      cout << "Residual[1] = " << Residual[1] << endl;
      LinSysRes.AddBlock(iPoint, Residual);

//      cout << "Implicit part -AA" << endl;
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }

// Diriclet BC
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Impose boundary values (Dirichlet) ---*/
//      Solution[0] = 0.0;
//      Solution[1] = 0.0;
//      node[iPoint]->SetSolution_Old(Solution);
//      LinSysRes.SetBlock_Zero(iPoint);
//
//      /*--- includes 1 in the diagonal ---*/
//      for (iVar = 0; iVar < nVar; iVar++) {
//        total_index = iPoint*nVar+iVar;
//        Jacobian.DeleteValsRowi(total_index);
//      }
//    }
//  }
  
  delete [] U_domain;
  delete [] U_wall;
  delete [] Normal;
  delete [] Residual;
  
}

void CTransLMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned long total_index;
  unsigned short iVar;

  //cout << "Arrived in BC_Far_Field. -AA" << endl;
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Impose boundary values (Dirichlet) ---*/
      Solution[0] = Intermittency_Inf;
      Solution[1] = REth_Inf;
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- includes 1 in the diagonal ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }

}

void CTransLMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker) {
//cout << "BC_Inlet reached -AA" << endl;
BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker) {
//cout << "BC_Outlet reached -AA" << endl;
BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker) {
//cout << "BC_Outlet reached -AA" << endl;
BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}
