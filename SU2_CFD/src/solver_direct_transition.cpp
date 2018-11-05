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

  /*--- Indicate that this is a solver for the transition model,
        not for the turbulence model. ---*/
  transitionSolver = true;
}

CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver() {

  /*--- Array initialization ---*/
  constants = NULL;

  /*--- Indicate that this is a solver for the transition model,
        not for the turbulence model. ---*/
  transitionSolver = true;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Dimension of the problem --> dependent on the turbulence model. ---*/
  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/
  nVarGrad = nVar;

  /*--- Determine the number of variables in the corresponding turbulence model. ---*/
  nVarTurbModel = 0;
  switch (config->GetKind_Turb_Model()) {
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      nVarTurbModel = 1;
      break;
    case SST:
      nVarTurbModel = 2;
      break;
  }

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
  constants = new su2double[11];
  constants[0] =  2.0;    // ca1
  constants[1] =  0.06;   // ca2
  constants[2]  =  1.0;    // ce1
  constants[3]  = 50.0;    // ce2
  constants[4]  =  0.03;   // cthetat
  constants[5]  =  2.0;    // s1
  constants[6]  =  1.0;    // sigmaf
  constants[7]  =  2.0;    // sigmathetat
  constants[8]  =  5.0;    // Flength_CF
  constants[9]  =  0.7;    // C_Fonset1_CF
  constants[10] =  0.6944; // CHe_max

  /*--- If the cross flow instability term must not be added, set the
        value of Flength_CF, constants[8], to zero. ---*/
  if (!config->GetLM_Cross_Flow_Instability()) constants[8] = 0.0;

  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];

  lowerlimit[0] = 0.0;
  upperlimit[0] = 2.0;

  lowerlimit[1] = 20.0;
  upperlimit[1] = 1.0e5;

  /*--- Store the far-field values of the intermittency and Re_theta. */
  Intermittency_Inf = 1.0;
  REth_Inf          = CSourcePieceWise_TransLM::GetREth(config->GetTurbulenceIntensity_FreeStream());

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

  /*--- Adjust the number of solution variables in the incompressible
   restart. We always carry a space in nVar for the energy equation in the
   mean flow solver, but we only write it to the restart if it is active.
   Therefore, we must reduce skipVars here if energy is inactive so that
   the turbulent variables are read correctly. ---*/
  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();

  if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

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

  unsigned long ExtIter = config->GetExtIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_turb     = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));

  /*--- Initialize the residual vector ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {
    LinSysRes.SetBlock_Zero(iPoint);
  }

  /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  if (limiter_turb) SetSolution_Limiter(geometry, config);

  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
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
    const su2double Re_theta_crit = CSourcePieceWise_TransLM::GetREth_crit(Re_theta);
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

  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Primitive variables w/o reconstruction and their gradients. ---*/
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

    /*--- Set vorticity and strain rate magnitude ---*/
    numerics->SetVorticity(solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(), NULL);
    numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(), 0.0);

    /*--- Laminar and eddy viscosity ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
    numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);

    /*--- Turbulence and transition variables. No gradients are needed for the source terms. ---*/
    numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), NULL);
    numerics->SetTransVar(node[iPoint]->GetSolution(), NULL);

    /*--- Set the volume and the distance to the wall. ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Compute the source term ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the Jacobian ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
  }
}

void CTransLMSolver::BC_Far_Field(CGeometry      *geometry,
                                  CSolver        **solver_container,
                                  CNumerics      *conv_numerics,
                                  CNumerics      *visc_numerics,
                                  CConfig        *config,
                                  unsigned short val_marker) {

  /*--- Easier storage whether or not grid movement takes place. ---*/
  const bool grid_movement = config->GetGrid_Movement();

  /*--- Loop over the vertices of this domain and check if
        it belongs to this domain (i.e. not a halo node). ---*/
  for (unsigned long iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const unsigned long iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Set the value of the flow solution at the infinity
            as well as at the node of the farfield boundary. ---*/
      su2double *V_infty  = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      su2double *V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Do the same for the transition variables. ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_i[1] = node[iPoint]->GetSolution(1);

      Solution_j[0] = Intermittency_Inf;
      Solution_j[1] = REth_Inf;

      /*--- Get Normal (it is necessary to change the sign) ---*/
      su2double Normal[] = {0.0, 0.0, 0.0};
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      Normal[0] = -Normal[0]; Normal[1] = -Normal[1]; Normal[2] = -Normal[2];

      /*--- Store the necessary info in conv_numerics, including the
            grid velocities, if necessary. Note that transition variables
            are stored as turbulent variables. ---*/
      conv_numerics->SetPrimitive(V_domain, V_infty);
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute residuals and Jacobians ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Add residuals and Jacobians ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
}

void CTransLMSolver::BC_Inlet(CGeometry      *geometry,
                              CSolver        **solver_container,
                              CNumerics      *conv_numerics,
                              CNumerics      *visc_numerics,
                              CConfig        *config,
                              unsigned short val_marker) {

  /*--- Call BC_Far_Field to do the actual work. ---*/
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_Outlet(CGeometry      *geometry,
                               CSolver        **solver_container,
                               CNumerics      *conv_numerics,
                               CNumerics      *visc_numerics,
                               CConfig        *config,
                               unsigned short val_marker) {

  /*--- Call BC_Far_Field to do the actual work. ---*/
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::BC_HeatFlux_Wall(CGeometry      *geometry,
                                      CSolver        **solver_container,
                                      CNumerics      *numerics,
                                      CNumerics      *visc_numerics,
                                      CConfig        *config,
                                      unsigned short val_marker) {

  /*--- Convective and diffusive fluxes are zero for the transition
        variables on a wall boundary. ---*/
}

void CTransLMSolver::BC_Isothermal_Wall(CGeometry      *geometry,
                                        CSolver        **solver_container,
                                        CNumerics      *conv_numerics,
                                        CNumerics      *visc_numerics,
                                        CConfig        *config,
                                        unsigned short val_marker) {

  /*--- Convective and diffusive fluxes are zero for the transition
        variables on a wall boundary. ---*/
}

void CTransLMSolver::BC_Inlet_MixingPlane(CGeometry      *geometry,
                                          CSolver        **solver_container,
                                          CNumerics      *conv_numerics,
                                          CNumerics      *visc_numerics,
                                          CConfig        *config,
                                          unsigned short val_marker) {

  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

void CTransLMSolver::BC_Inlet_Turbo(CGeometry      *geometry,
                                    CSolver        **solver_container,
                                    CNumerics      *conv_numerics,
                                    CNumerics      *visc_numerics,
                                    CConfig        *config,
                                    unsigned short val_marker) {

  /*--- Get the required information from config. ---*/
  const unsigned short nSpanWiseSections = config->GetnSpanWiseSections();
  const bool           grid_movement     = config->GetGrid_Movement();

  /*--- Loop over the number of spanwise sections. ---*/
  for (unsigned short iSpan= 0; iSpan < nSpanWiseSections ; iSpan++) {

    /*--- Loop over all the vertices on this boundary marker ---*/
    for (long iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- Find the node related to the vertex ---*/
      const unsigned long iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const unsigned long oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const unsigned long Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      su2double Normal[3];
      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Retrieve the primitive variables and set them in the numerics class. ---*/
      su2double *V_inlet  = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);
      su2double *V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Store the transition variables in the turbulence variables. The reason
            is that a lot of work of the transition solver is carried out by its
            base class, which is the base class of the turbulence solvers. ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_i[1] = node[iPoint]->GetSolution(1);

      Solution_j[0] = Intermittency_Inf;
      Solution_j[1] = REth_Inf;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set the normal and the grid velocities, if needed. ---*/
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Primitive variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    }
  }
}

void CTransLMSolver::BC_Fluid_Interface(CGeometry *geometry,
                                        CSolver   **solver_container,
                                        CNumerics *conv_numerics,
                                        CNumerics *visc_numerics,
                                        CConfig   *config) {

  /*--- Determine whether or not grid movement is present and
        the number of primitive flow variables. ---*/
  const bool grid_movement      = config->GetGrid_Movement();
  const unsigned short nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();

  /*--- Allocate the memory of some local variables. ---*/
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];

  /*--- Loop over the markers and select the fluid interfaces. ---*/
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      /*--- Loop over the vertices of this marker. ---*/
      for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Get the point ID and determine if it is owned by this rank. ---*/
        const unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Determine the number of donor vertices for this point. ---*/
          const unsigned long nDonorVertex = GetnSlidingStates(iMarker, iVertex);
 
          /*--- Initialize Residual, this will serve to accumulate the average ---*/
          for (unsigned short iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertices and compute the averaged flux ---*/
          su2double Normal[3];
          for (unsigned long jVertex = 0; jVertex < nDonorVertex; jVertex++){

            /*--- Get the normal and negate it. ---*/
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

            /*--- Get the primitive flow variables. ---*/
            for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }

            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/
            const su2double weight = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set the primitive and transition (turbulent) variables. ---*/
            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );

            Solution_i[0] = node[iPoint]->GetSolution(0);
            Solution_i[1] = node[iPoint]->GetSolution(1);

            Solution_j[0] = GetSlidingState(iMarker, iVertex, 0, jVertex);
            Solution_j[1] = GetSlidingState(iMarker, iVertex, 1, jVertex);

            conv_numerics->SetTurbVar(Solution_i, Solution_j);

            /*--- Set the normal vector and the grid velocities, if needed. ---*/
            conv_numerics->SetNormal(Normal);

            if (grid_movement)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

            /*--- Compute the residual, which is stored in tmp_residual. ---*/
            conv_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);

            /*--- Accumulate the residuals to compute the average ---*/
            for (unsigned short iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] += weight*tmp_residual[iVar];
          }

          /*--- Add Residuals and Jacobians ---*/
          LinSysRes.AddBlock(iPoint, Residual);
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          /*--- Get the point ID for the coordinates of the normal vector. ---*/
          const unsigned long Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          /*--- Set the normal vector and the coordinates ---*/
          visc_numerics->SetNormal(Normal);
          visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

          /*--- Primitive variables ---*/
          visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);

          /*--- Transition (turbulent) variables and its gradients  ---*/
          visc_numerics->SetTurbVar(Solution_i, Solution_j);
          visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

          /*--- Compute and update residual ---*/
          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] tmp_residual;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
}


void CTransLMSolver::SetInletAtVertex(su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {

  const unsigned short intermit_pos = nDim+2+nDim+nVarTurbModel;
  const unsigned short Re_theta_pos = intermit_pos+1;

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[intermit_pos];
  Inlet_TurbVars[iMarker][iVertex][1] = val_inlet[Re_theta_pos];
}

su2double CTransLMSolver::GetInletAtVertex(su2double *val_inlet,
                                           unsigned long val_inlet_point,
                                           unsigned short val_kind_marker,
                                           string val_marker,
                                           CGeometry *geometry,
                                           CConfig *config) {

  /*--- Initialize the return value to zero. ---*/
  su2double Area = 0.0;

  /*--- Test for an inlet. ---*/
  if (val_kind_marker == INLET_FLOW) {

    /*--- Alias positions within inlet file for readability ---*/
    const unsigned short intermit_pos = nDim+2+nDim+nVarTurbModel;
    const unsigned short Re_theta_pos = intermit_pos+1;

    /*--- Loop over all markers and select the correct one. ---*/
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {

        /*--- Loop over all the vertices belonging to this marker. ---*/
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

          /*--- Check for the correct inlet point. ---*/
          const unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (iPoint == val_inlet_point) {

            /*-- Compute boundary face area for this vertex. ---*/
            su2double Normal[3];
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = 0.0;
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);

            /*--- Access and store the inlet variables for this vertex. ---*/
            val_inlet[intermit_pos]  = Inlet_TurbVars[iMarker][iVertex][0];
            val_inlet[Re_theta_pos] = Inlet_TurbVars[iMarker][iVertex][1];

            /*--- Exit once we find the point. ---*/
            return Area;
          }
        }
      }
    }
  }

  /*--- If we don't find a match, then the child point is not on the
        current inlet boundary marker. Return zero area so this point does
        not contribute to the restriction operator and continue. ---*/
  return Area;
}

void CTransLMSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {

  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_TurbVars[iMarker][iVertex][0] = Intermittency_Inf;
    Inlet_TurbVars[iMarker][iVertex][1] = REth_Inf;
  }
}
