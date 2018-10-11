/*!
 * \file solution_direct_transition.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author A. Aranake
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

CTransLMSolver::CTransLMSolver(void) : CTurbSolver() {

  /*--- Array initialization ---*/
  constants = NULL;
  Inlet_TurbVars = NULL;
}

CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver() {

  /*--- Array initialization ---*/
  constants = NULL;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Dimension of the problem --> dependent on the turbulence model. ---*/
  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/
  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];

  /*--- Single grid simulation ---*/ 
  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/
    Residual = new su2double[nVar];     for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
    Residual_RMS = new su2double[nVar]; for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_i = new su2double[nVar];   for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
    Residual_j = new su2double[nVar];   for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
    Residual_Max = new su2double[nVar]; for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;

    /*--- Define some structures for locating max residuals ---*/
    Point_Max = new unsigned long[nVar];
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (unsigned short iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }

    /*--- Define some auxiliary vector related with the solution ---*/
    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

    /*--- Define some auxiliary vector related with the geometry ---*/
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];

    /*--- Define some auxiliary vector related with the flow solution ---*/
    FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];

    /*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (LM model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      const unsigned short nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }

  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    /*--- S matrix := inv(R)*transpose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];

    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }

  /*--- Initialize the values for the model constants ---*/
  constants = new su2double[8];
  constants[0] =  2.0;    // ca1
  constants[1] =  0.06;   // ca2
  constants[2] =  1.0;    // ce1
  constants[3] = 50.0;    // ce2
  constants[4] =  0.03;   // cthetat
  constants[5] =  2.0;    // s1
  constants[6] =  1.0;    // sigmaf
  constants[7] =  2.0;    // sigmathetat

  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];

  lowerlimit[0] = 0.0;
  upperlimit[0] = 2.0;

  lowerlimit[1] = 20.0;
  upperlimit[1] = 1.0e5;

  /*--- Store the far-field values of the intermittency and Re_theta. */
  Intermittency_Inf = config->GetIntermittency_FreeStream();
  REth_Inf          = GetREth(config->GetTurbulenceIntensity_FreeStream());

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CTransLMVariable(Intermittency_Inf, REth_Inf, nDim, nVar, config);

  /*--- MPI solution ---*/

//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);

  /*--- Initializate quantities for SlidingMesh Interface ---*/
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {

    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (unsigned long iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++) {
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (unsigned short iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
       due to arbitrary number of turbulence variables ---*/
  Inlet_TurbVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker] = new su2double*[nVertex[iMarker]];
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_TurbVars[iMarker][iVertex] = new su2double[nVar];
      Inlet_TurbVars[iMarker][iVertex][0] = Intermittency_Inf;
      Inlet_TurbVars[iMarker][iVertex][1] = REth_Inf;
    }
  }
}

CTransLMSolver::~CTransLMSolver(void) {
 
  if (constants != NULL) delete [] constants;

  if (SlidingState != NULL) {
    for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      if (SlidingState[iMarker] != NULL) {
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if (SlidingState[iMarker][iVertex] != NULL) {
            for (unsigned short iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }

  if (SlidingStateNodes != NULL) {
    for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      if (SlidingStateNodes[iMarker] != NULL)
        delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }
}

su2double CTransLMSolver::GetREth_crit(const su2double var_Re_theta) {

  su2double Recrit;

  /* Menter correlation, also mentioned at https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html. */
  if(var_Re_theta <= 1870.0) {

   /* Quartic expression is used. Define the constants. */
   const su2double a0 = -3.96035;
   const su2double a1 =  1.0120656;
   const su2double a2 = -8.68230e-4;
   const su2double a3 =  6.96506e-7;
   const su2double a4 = -1.74105e-10;

   /* Compute the value of Recrit. */
   const su2double val1 = var_Re_theta;
   const su2double val2 = val1*val1;
   const su2double val3 = val1*val2;
   const su2double val4 = val2*val2;

   Recrit = a0 + a1*val1 + a2*val2 + a3*val3 + a4*val4;

  } else {

    /* Use correlation valid for Re_theta larger than 1870.0. */
    Recrit = var_Re_theta - (593.11 + 0.482*(var_Re_theta-1870.0));
  }

  /* Malan correlation. */
  // Recrit = min(0.615*var_Re_theta+61.5, var_Re_theta);

  /* Return the value of Recrit. */
  return Recrit;
}

void CTransLMSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  string restart_filename = config->GetSolution_FlowFileName();

  /*--- Modify file name for multizone problems ---*/
  if (config->GetnZone() > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, config->GetiZone());

  /*--- Modify file name for an unsteady restart ---*/
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);

  if (dual_time|| time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
  if (config->GetRead_Binary_Restart())
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  else
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);

  /*--- Skip flow and turbulence variables. ---*/
  unsigned short skipVars = solver[MESH_0][FLOW_SOL]->GetnVar()
                          + solver[MESH_0][TURB_SOL]->GetnVar() + nDim;

  /*--- Load data from the restart into correct containers. ---*/
  unsigned long iPoint_Global_Local = 0;
  for (unsigned long iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    long iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/
      unsigned long index = iPoint_Global_Local*Restart_Vars[1] + skipVars;
      for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  if (iPoint_Global_Local < nPointDomain) sbuf_NotMatching = 1;

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);

  /*--- MPI solution and compute the separation intermittency. ---*/

//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart.
  solver[MESH_0][TRANS_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][TRANS_SOL]->Set_MPI_Solution(geometry[MESH_0], config);

  solver[MESH_0][TRANS_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

  /*--- I don't think it is necessary to interpolate the transition variables to the coarser
        grids, because the transition and turbulence equations are only solved on the finest
        grid and the transition variables only interact with the turbulence solver, not with
        the mean flow solver. ---*/
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

  /*--- Determine whether a Spalart-Allmaras type turbulence model is used. ---*/
  const unsigned short turbModel = config->GetKind_Turb_Model();
  const bool SA_turb  = (turbModel == SA) || (turbModel == SA_NEG) || (turbModel == SA_E) ||
                        (turbModel == SA_COMP) || (turbModel == SA_E_COMP);

  /*--- Easier storage of some of the model constants. */
  const su2double ce2 = constants[3];
  const su2double s1  = constants[5];

  /*--- Compute the value of gamma_sep, which models separation-induced transition. ---*/
  for (unsigned int iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {

    /*--- Get the required variables for the mean flow, turbulence and transition
          variables. Note that it is not needed to recompute the gradients of the
          flow variables, because these have not been altered (assuming that the
          transition solver is called after the turbulence solver). ---*/
    const su2double dist = geometry->node[iPoint]->GetWall_Distance();

    const su2double rho        = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    const su2double vel2       = solver_container[FLOW_SOL]->node[iPoint]->GetVelocity2();
    const su2double mu         = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    const su2double strMag     = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    const su2double *vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    const su2double Omega      = sqrt(vorticity[0]*vorticity[0] + vorticity[1]*vorticity[1]
                               +      vorticity[2]*vorticity[2]);

    const su2double muT   = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
    const su2double kine  = SA_turb ? 0.0 : solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    const su2double omega = SA_turb ? 1.0 : solver_container[TURB_SOL]->node[iPoint]->GetSolution(1);

    const su2double intermittency = node[iPoint]->GetSolution(0);
    const su2double Re_theta      = node[iPoint]->GetSolution(1);

    /*--- Compute the different Reynolds numbers that appear in the formulation. ---*/
    const su2double Rev           = rho*strMag*dist*dist/mu;
    const su2double Re_theta_crit = GetREth_crit(Re_theta);
    const su2double Re_omega      = rho*omega*dist*dist/mu;

    /*--- Compute the values of Fwake and RT. Here a distinction must be made between
          Spalart-Allmaras type models and the SST model. ---*/
    su2double Fwake, RT;
    if( SA_turb ) {
      Fwake = 1.0;
      RT    = muT/mu;
    } else {
      const su2double val = 1.e-5*Re_omega;
      Fwake = exp(-val*val);
      RT    = rho*kine/(mu*omega);
    }

    /*--- Compute the ratio d over delta and Freattach. ---*/
    const su2double dOverDelta = rho*vel2/(375.0*Omega*mu*Re_theta);
    const su2double valR       = 0.05*RT;
    const su2double FReattach  = exp(-valR*valR*valR*valR);

    /*--- Compute the value of Ftheta_t. */
    su2double val1 = Fwake*exp(-dOverDelta*dOverDelta*dOverDelta*dOverDelta);
    su2double val2 = (ce2*intermittency-1.0)/(ce2-1.0);
    su2double val3 = 1.0 - val2*val2;

    const su2double Ftheta_t = min(max(val1,val3),1.0);

    /*--- Compute the value of gamma_sep and store it. ---*/
    val1 = Rev/(3.235*Re_theta_crit) - 1.0;
    val2 = max(0.0,val1)*FReattach;
    val3 = min(s1*val2,2.0)*Ftheta_t;

    node[iPoint]->SetGammaSep(val3);
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

  delete [] U_domain;
  delete [] U_wall;
  delete [] Normal;
  delete [] Residual;
  
}

void CTransLMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  unsigned long total_index;
  unsigned short iVar;

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
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker) {
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker) {
  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker); 
}
