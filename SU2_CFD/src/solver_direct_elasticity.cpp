/*!
 * \file solver_direct_elasticity.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
 * \author R. Sanchez
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

CFEM_ElasticitySolver::CFEM_ElasticitySolver(void) : CSolver() {
  
  nElement = 0;
  nDim = 0;
  nMarker = 0;
  
  nFEA_Terms = 1;
  
  nPoint = 0;
  nPointDomain = 0;
  
  Total_CFEA = 0.0;
  WAitken_Dyn = 0.0;
  WAitken_Dyn_tn1 = 0.0;
  loadIncrement = 1.0;
  
  element_container = NULL;
  node = NULL;
  
  element_properties = NULL;
  elProperties = NULL;
  
  GradN_X = NULL;
  GradN_x = NULL;
  
  Jacobian_c_ij = NULL;
  Jacobian_s_ij = NULL;
  Jacobian_k_ij = NULL;
  
  MassMatrix_ij = NULL;
  
  mZeros_Aux = NULL;
  mId_Aux = NULL;
  
  Res_Stress_i = NULL;
  Res_Ext_Surf = NULL;
  Res_Time_Cont = NULL;
  Res_FSI_Cont = NULL;
  
  Res_Dead_Load = NULL;
  
  nodeReactions = NULL;
  
  solutionPredictor = NULL;
  
  SolRest = NULL;
  
  normalVertex = NULL;
  stressTensor = NULL;
  
  Solution_Interm = NULL;
  
  iElem_iDe   = NULL;
  
}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
  unsigned long iPoint;
  unsigned short iVar, jVar, iDim, jDim;
  unsigned short iTerm, iKind;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation
  bool gen_alpha = (config->GetKind_TimeIntScheme_FEA() == GENERALIZED_ALPHA);  // Generalized alpha method requires residual at previous time step.
  
  bool de_effects = config->GetDE_Effects();                      // Test whether we consider dielectric elastomers
  
  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
  bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
  
  element_based = false;          // A priori we don't have an element-based input file (most of the applications will be like this)
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();
  
  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();
  
  /*--- Number of different terms for FEA ---*/
  nFEA_Terms = 1;
  if (de_effects) nFEA_Terms++;       // The DE term is DE_TERM = 1
  if (incompressible) nFEA_Terms = 3; // The incompressible term is INC_TERM = 2
  
  /*--- Here is where we assign the kind of each element ---*/
  
  /*--- First level: different possible terms of the equations ---*/
  element_container = new CElement** [MAX_TERMS];
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS];
  
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++) {
    for (iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
      element_container[iTerm][iKind] = NULL;
    }
  }
  
  if (nDim == 2) {
    
      /*--- Basic terms ---*/
      element_container[FEA_TERM][EL_TRIA] = new CTRIA1(nDim, config);
      element_container[FEA_TERM][EL_QUAD] = new CQUAD4(nDim, config);
    
      if (de_effects){
        element_container[DE_TERM][EL_TRIA] = new CTRIA1(nDim, config);
        element_container[DE_TERM][EL_QUAD] = new CQUAD4(nDim, config);
      }

      if (incompressible){
        element_container[INC_TERM][EL_TRIA] = new CTRIA1(nDim, config);
        element_container[INC_TERM][EL_QUAD] = new CQUAD1(nDim, config);
      }

  }
  else if (nDim == 3) {

      element_container[FEA_TERM][EL_TETRA] = new CTETRA1(nDim, config);
      element_container[FEA_TERM][EL_HEXA] = new CHEXA8(nDim, config);

      if (de_effects){
        element_container[DE_TERM][EL_TETRA] = new CTETRA1(nDim, config);
        element_container[DE_TERM][EL_HEXA] = new CHEXA8(nDim, config);
      }

      if (incompressible) {
        element_container[INC_TERM][EL_TETRA] = new CTETRA1(nDim, config);
        element_container[INC_TERM][EL_HEXA] = new CHEXA1(nDim, config);
      }


  }
  
  node              = new CVariable*[nPoint];
  
  /*--- Set element properties ---*/
  elProperties = new unsigned long[4];
  for (iVar = 0; iVar < 4; iVar++)
    elProperties[iVar] = 0;
  Set_ElementProperties(geometry, config);
  
  GradN_X = new su2double [nDim];
  GradN_x = new su2double [nDim];
  
  Total_CFEA      = 0.0;
  WAitken_Dyn       = 0.0;
  WAitken_Dyn_tn1   = 0.0;
  loadIncrement     = 0.0;
  
  SetFSI_ConvValue(0,0.0);
  SetFSI_ConvValue(1,0.0);
  
  nVar = nDim;
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual = new su2double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Residual i and residual j for the Matrix-Vector Product for the Mass Residual ---*/
  if (dynamic){
    Residual_i = new su2double[nVar];        for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
    Residual_j = new su2double[nVar];        for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
  
  Solution_Interm = NULL;
  if (gen_alpha) {
    Solution_Interm = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Solution_Interm[iVar] = 0.0;
  }
  
  nodeReactions = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) nodeReactions[iVar]   = 0.0;
  
  /*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/
  
  unsigned short nSolVar;
  string text_line, filename;
  ifstream restart_file;

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;
  
  SolRest = new su2double[nSolVar];

  /*--- Initialize from zero everywhere. ---*/

  for (iVar = 0; iVar < nSolVar; iVar++) SolRest[iVar] = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
  }
  
  bool reference_geometry = config->GetRefGeom();
  if (reference_geometry) Set_ReferenceGeometry(geometry, config);
  
  bool prestretch_fem = config->GetPrestretch();
  if (prestretch_fem) Set_Prestretch(geometry, config);  
  
  /*--- Term ij of the Jacobian ---*/
  
  Jacobian_ij = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_ij[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_ij[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Term ij of the Mass Matrix (only if dynamic analysis) ---*/
  MassMatrix_ij = NULL;
  if (dynamic) {
    MassMatrix_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      MassMatrix_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        MassMatrix_ij[iVar][jVar] = 0.0;
      }
    }
  }
  
  Jacobian_c_ij = NULL;
  Jacobian_s_ij = NULL;
  if (nonlinear_analysis) {
    
    /*--- Term ij of the Jacobian (constitutive contribution) ---*/
    
    Jacobian_c_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_c_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_c_ij[iVar][jVar] = 0.0;
      }
    }
    
    /*--- Term ij of the Jacobian (stress contribution) ---*/
    
    Jacobian_s_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_s_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_s_ij[iVar][jVar] = 0.0;
      }
    }
    
  }
  
  /*--- Term ij of the Jacobian (incompressibility term) ---*/
  Jacobian_k_ij = NULL;
  if (incompressible) {
    Jacobian_k_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_k_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_k_ij[iVar][jVar] = 0.0;
      }
    }
  }
  
  /*--- Stress contribution to the node i ---*/
  Res_Stress_i = new su2double[nVar];
  
  /*--- Contribution of the external surface forces to the residual (auxiliary vector) ---*/
  Res_Ext_Surf = new su2double[nVar];
  
  /*--- Contribution of the body forces to the residual (auxiliary vector) ---*/
  Res_Dead_Load = NULL;
  if (body_forces) {
    Res_Dead_Load = new su2double[nVar];
  }
  
  /*--- Contribution of the fluid tractions to the residual (auxiliary vector) ---*/
  Res_FSI_Cont = NULL;
  if (fsi) {
    Res_FSI_Cont = new su2double[nVar];
  }
  
  /*--- Time integration contribution to the residual ---*/
  Res_Time_Cont = NULL;
  if (dynamic) {
    Res_Time_Cont = new su2double [nVar];
  }
  
  /*--- Matrices to impose clamped boundary conditions ---*/
  
  mZeros_Aux = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++)
    mZeros_Aux[iDim] = new su2double[nDim];
  
  mId_Aux = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++)
    mId_Aux[iDim] = new su2double[nDim];
  
  for(iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      mZeros_Aux[iDim][jDim] = 0.0;
      mId_Aux[iDim][jDim] = 0.0;
    }
    mId_Aux[iDim][iDim] = 1.0;
  }
  
  
  /*--- Initialization of matrix structures ---*/
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;
  
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  
  if (dynamic) {
    MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
    TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
    TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
  /*--- Initialization of linear solver structures ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  LinSysReact.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Initialize the auxiliary vector and matrix for the computation of the nodal Reactions ---*/
  
  normalVertex = new su2double [nDim];
  
  stressTensor = new su2double* [nDim];
  for (iVar = 0; iVar < nVar; iVar++) {
    stressTensor[iVar] = new su2double [nDim];
  }
  
  /*---- Initialize the auxiliary vector for the solution predictor ---*/
  
  solutionPredictor = new su2double [nVar];
  
  iElem_iDe = NULL;

  /*--- Initialize the value of the total objective function ---*/
   Total_OFRefGeom = 0.0;
   Total_OFRefNode = 0.0;

   /*--- Initialize the value of the global objective function ---*/
   Global_OFRefGeom = 0.0;
   Global_OFRefNode = 0.0;

   /*--- Initialize the value of the total gradient for the forward mode ---*/
   Total_ForwardGradient = 0.0;

   if (config->GetDirectDiff() == D_YOUNG ||
       config->GetDirectDiff() == D_POISSON ||
       config->GetDirectDiff() == D_RHO ||
       config->GetDirectDiff() == D_RHO_DL ||
       config->GetDirectDiff() == D_EFIELD ||
       config->GetDirectDiff() == D_MACH ||
       config->GetDirectDiff() == D_PRESSURE){

     /*--- Header of the temporary output file ---*/
     ofstream myfile_res;

     if (config->GetDirectDiff() == D_YOUNG) myfile_res.open ("Output_Direct_Diff_E.txt");
     if (config->GetDirectDiff() == D_POISSON) myfile_res.open ("Output_Direct_Diff_Nu.txt");
     if (config->GetDirectDiff() == D_RHO) myfile_res.open ("Output_Direct_Diff_Rho.txt");
     if (config->GetDirectDiff() == D_RHO_DL) myfile_res.open ("Output_Direct_Diff_Rho_DL.txt");
     if (config->GetDirectDiff() == D_EFIELD) myfile_res.open ("Output_Direct_Diff_EField.txt");

     if (config->GetDirectDiff() == D_MACH) myfile_res.open ("Output_Direct_Diff_Mach.txt");
     if (config->GetDirectDiff() == D_PRESSURE) myfile_res.open ("Output_Direct_Diff_Pressure.txt");

     myfile_res << "Objective Function " << "\t";

     if (dynamic) myfile_res << "O. Function Averaged " << "\t";

     myfile_res << "Sensitivity Local" << "\t";

     myfile_res << "Sensitivity Averaged" << "\t";

     myfile_res << endl;

     myfile_res.close();

   }

   /*--- Penalty value - to maintain constant the stiffness in optimization problems - TODO: this has to be improved ---*/
   PenaltyValue = 0.0;

  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- If dynamic, we also need to communicate the old solution ---*/
  
  if(dynamic) Set_MPI_Solution_Old(geometry, config);
  
}

CFEM_ElasticitySolver::~CFEM_ElasticitySolver(void) {
  
  unsigned short iVar, jVar;
  unsigned long iElem;
  
  if (element_container != NULL) {
    for (iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        if (element_container[iVar][jVar] != NULL) delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }
  
  if (element_properties != NULL){
    for (iElem = 0; iElem < nElement; iElem++)
      if (element_properties[iElem] != NULL) delete element_properties[iElem];
    delete [] element_properties;
  }
  
  for (iVar = 0; iVar < nVar; iVar++) {
    if (Jacobian_s_ij != NULL) delete [] Jacobian_s_ij[iVar];
    if (Jacobian_c_ij != NULL) delete [] Jacobian_c_ij[iVar];
    if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij[iVar];
    delete [] mZeros_Aux[iVar];
    delete [] mId_Aux[iVar];
    delete [] stressTensor[iVar];
  }
  
  if (Jacobian_s_ij != NULL) delete [] Jacobian_s_ij;
  if (Jacobian_c_ij != NULL) delete [] Jacobian_c_ij;
  if (Jacobian_k_ij != NULL) delete [] Jacobian_k_ij;
  delete [] Res_Stress_i;
  delete [] Res_Ext_Surf;
  if (Res_Time_Cont != NULL) delete[] Res_Time_Cont;
  if (Res_Dead_Load != NULL) delete[] Res_Dead_Load;
  delete [] SolRest;
  delete [] GradN_X;
  delete [] GradN_x;
  
  delete [] mZeros_Aux;
  delete [] mId_Aux;
  
  delete [] nodeReactions;
  
  delete [] normalVertex;
  delete [] stressTensor;
  
  if (iElem_iDe != NULL) delete [] iElem_iDe;
}

void CFEM_ElasticitySolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  
  unsigned short nSolVar;
  
  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
        if (dynamic) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
            Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
          }
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        if (dynamic) {
          for (iVar = nVar; iVar < 3*nVar; iVar++)
            Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        }
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);
        
        if (dynamic) {
          
          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->SetSolution_Vel(iVar, SolRest[iVar+nVar]);
            node[iPoint]->SetSolution_Accel(iVar, SolRest[iVar+2*nVar]);
          }
          
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  unsigned short nSolVar;
  
  nSolVar = 3 * nVar;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_time_n(iVar);
          Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel_time_n(iVar);
          Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel_time_n(iVar);
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nSolVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          node[iPoint]->SetSolution_time_n(iVar, SolRest[iVar]);
          node[iPoint]->SetSolution_Vel_time_n(iVar, SolRest[iVar+nVar]);
          node[iPoint]->SetSolution_Accel_time_n(iVar, SolRest[iVar+2*nVar]);
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_DispOnly(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;         nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_Pred(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;     nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy predicted solution variables back into the variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_Pred_Old(CGeometry *geometry, CConfig *config) {
  
  /*--- We are communicating the solution predicted, current and old, and the old solution ---*/
  /*--- necessary for the Aitken relaxation ---*/
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  /*--- Analogous to the dynamic solution, in this case we need 3 * nVar variables per node ---*/
  unsigned short nSolVar;
  nSolVar = 3 * nVar;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
          Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
          Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred_Old(iVar);
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nSolVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[(iVar+nVar)*nVertexR+iVertex]);
          node[iPoint]->SetSolution_Pred_Old(iVar, Buffer_Receive_U[(iVar+2*nVar)*nVertexR+iVertex]);
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) {

  unsigned long iElem;
  unsigned long index;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  string filename;
  ifstream properties_file;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  element_properties = new CElementProperty*[nElement];

  /*--- Restart the solution from file information ---*/

  filename = config->GetFEA_FileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no element-based properties file. All elements get uniform properties." << endl;

    for (iElem = 0; iElem < nElement; iElem++){
      element_properties[iElem] = new CElementProperty(0, 0, 0, 0);
    }

    element_based = false;

  }
  else{

    element_based = true;

    /*--- In case this is a parallel simulation, we need to perform the
       Global2Local index transformation first. ---*/

    long *Global2Local = new long[geometry->GetGlobal_nElemDomain()];

    /*--- First, set all indices to a negative value by default ---*/

    for (iElem = 0; iElem < geometry->GetGlobal_nElemDomain(); iElem++)
      Global2Local[iElem] = -1;

    /*--- Now fill array with the transform values only for the points in the rank (including halos) ---*/

    for (iElem = 0; iElem < nElement; iElem++)
      Global2Local[geometry->elem[iElem]->GetGlobalIndex()] = iElem;

    /*--- Read all lines in the restart file ---*/

    long iElem_Local;
    unsigned long iElem_Global_Local = 0, iElem_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/

    getline (properties_file, text_line);

    for (iElem_Global = 0; iElem_Global < geometry->GetGlobal_nElemDomain(); iElem_Global++ ) {

      getline (properties_file, text_line);

      istringstream point_line(text_line);

      /*--- Retrieve local index. If this element from the restart file lives
         only on a different processor, the value of iPoint_Local will be -1.
         Otherwise, the local index for this node on the current processor
         will be returned and used to instantiate the vars. ---*/

      iElem_Local = Global2Local[iElem_Global];

      if (iElem_Local >= 0) {

        if (config->GetDE_Effects())
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];
        else
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];

        element_properties[iElem_Local] = new CElementProperty(elProperties[0],
                                                         elProperties[1],
                                                         elProperties[2],
                                                         elProperties[3]);

        iElem_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iElem_Global_Local < nElement) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The properties file " << filename.data() << " doesn't match with the mesh file!" << endl;
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

    /*--- Close the restart file ---*/

    properties_file.close();

    /*--- Free memory needed for the transformation ---*/

    delete [] Global2Local;


  }


}


void CFEM_ElasticitySolver::Set_Prestretch(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint;
  unsigned long index;
  
  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  
  string filename;
  ifstream prestretch_file;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  
  /*--- Restart the solution from file information ---*/
  
  filename = config->GetPrestretch_FEMFileName();
  
  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);
  
  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;
  
  prestretch_file.open(filename.data(), ios::in);
  
  /*--- In case there is no file ---*/
  
  if (prestretch_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no FEM prestretch reference file!!" << endl;
    exit(EXIT_FAILURE);
  }
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  
  map<unsigned long,unsigned long> Global2Local;
  map<unsigned long,unsigned long>::const_iterator MI;
  
  /*--- Now fill array with the transform values only for local points ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
  
  /*--- Read all lines in the restart file ---*/
  
  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  
  /*--- The first line is the header ---*/
  
  getline (prestretch_file, text_line);
  
  while (getline (prestretch_file, text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    MI = Global2Local.find(iPoint_Global);
    if (MI != Global2Local.end()) {
    
      iPoint_Local = Global2Local[iPoint_Global];
      
      if (nDim == 2) point_line >> Solution[0] >> Solution[1] >> index;
      if (nDim == 3) point_line >> Solution[0] >> Solution[1] >> Solution[2] >> index;
      
      for (iVar = 0; iVar < nVar; iVar++) node[iPoint_Local]->SetPrestretch(iVar, Solution[iVar]);
      
      iPoint_Global_Local++;
    }
    iPoint_Global++;
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
  
  /*--- Close the restart file ---*/
  
  prestretch_file.close();
  
  /*--- We need to communicate here the prestretched geometry for the halo nodes. ---*/
  /*--- We avoid creating a new function as this may be reformatted.              ---*/

  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetPrestretch(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetPrestretch(iVar, SolRest[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }


}

void CFEM_ElasticitySolver::Set_ReferenceGeometry(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned long index;

  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  unsigned short file_format = config->GetRefGeom_FileFormat();

  string filename;
  su2double dull_val;
  ifstream reference_file;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  /*--- Restart the solution from file information ---*/

  filename = config->GetRefGeom_FEMFileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  reference_file.open(filename.data(), ios::in);

  /*--- In case there is no file ---*/

  if (reference_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no FEM reference geometry file!!" << endl;
    exit(EXIT_FAILURE);
  }

  if (rank == MASTER_NODE) cout << "Filename: " << filename << " and format " << file_format << "." << endl;

  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/

  long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];

  /*--- First, set all indices to a negative value by default ---*/

  for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
    Global2Local[iPoint] = -1;

  /*--- Now fill array with the transform values only for local points ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

  /*--- Read all lines in the restart file ---*/

  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- The first line is the header ---*/

  getline (reference_file, text_line);

  while (getline (reference_file, text_line)) {
    istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/

    iPoint_Local = Global2Local[iPoint_Global];

    if (iPoint_Local >= 0) {

      if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1];
      if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];

      for (iVar = 0; iVar < nVar; iVar++) node[iPoint_Local]->SetReference_Geometry(iVar, Solution[iVar]);

      iPoint_Global_Local++;
    }
    iPoint_Global++;
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

  /*--- I don't think we need to communicate ---*/

  /*--- Close the restart file ---*/

  reference_file.close();

  /*--- Free memory needed for the transformation ---*/

  delete [] Global2Local;

}



void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {
  
  
  unsigned long iPoint;
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);                          // Checks if it is the first iteration
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool restart = config->GetRestart();                        // Restart analysis
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter()); // Initial calculation for restart
  
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);     // Discrete adjoint FEM solver
  
  bool incremental_load = config->GetIncrementalLoad();               // If an incremental load is applied
  
  bool body_forces = config->GetDeadLoad();                     // Body forces (dead loads).

  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysAux.SetBlock_Zero(iPoint);
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }
  
  /*--- Set matrix entries to zero ---*/
  
  /*
   * If the problem is linear, we only need one Jacobian matrix in the problem, because
   * it is going to be constant along the calculations. Therefore, we only initialize
   * the Jacobian matrix once, at the beginning of the simulation.
   *
   * We don't need first_iter, because there is only one iteration per time step in linear analysis.
   */
  if ((initial_calc && linear_analysis)||
      (restart && initial_calc_restart && linear_analysis) ||
      (dynamic && disc_adj_fem)) {
    Jacobian.SetValZero();
  }
  
  /*
   * If the problem is dynamic, we need a mass matrix, which will be constant along the calculation
   * both for linear and nonlinear analysis. Only initialized once, at the first time step.
   *
   * The same with the integration constants, as for now we consider the time step to be constant.
   *
   * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
   */
  if ((dynamic && initial_calc && first_iter) ||
      (dynamic && restart && initial_calc_restart && first_iter) ||
      (dynamic && disc_adj_fem)) {
    MassMatrix.SetValZero();
    Compute_IntegrationConstants(config);
    Compute_MassMatrix(geometry, solver_container, numerics, config);
  }
  
  /*
   * If body forces are taken into account, we need to compute the term that goes into the residual,
   * which will be constant along the calculation both for linear and nonlinear analysis.
   *
   * Only initialized once, at the first iteration or the beginning of the calculation after a restart.
   *
   * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
   */
  
  if ((body_forces && initial_calc && first_iter) ||
      (body_forces && restart && initial_calc_restart && first_iter)) {
    // If the load is incremental, we have to reset the variable to avoid adding up over the increments
    if (incremental_load) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_BodyForces_Res();
    }
    // Compute the dead load term
    Compute_DeadLoad(geometry, solver_container, numerics, config);
  }
  
  /*
   * If the problem is nonlinear, we need to initialize the Jacobian and the stiffness matrix at least at the beginning
   * of each time step. If the solution method is Newton Rapshon, we initialize it also at the beginning of each
   * iteration.
   */
  
  if ((nonlinear_analysis) && ((newton_raphson) || (first_iter)))  {
    Jacobian.SetValZero();
    //    StiffMatrix.SetValZero();
  }
  
  /*
   * Some external forces may be considered constant over the time step.
   */
  if (first_iter)  {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_SurfaceLoad_Res();
  }
  
  /*
   * If we apply nonlinear forces, we need to clear the residual on each iteration
   */
  unsigned short iMarker;
  unsigned long iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case LOAD_BOUNDARY:
        /*--- Only if the load is nonzero - reduces computational cost ---*/
        if(config->GetLoad_Value(config->GetMarker_All_TagBound(iMarker)) != 0 ) {
          /*--- For all the vertices in the marker iMarker ---*/
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            /*--- Retrieve the point ID ---*/
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            /*--- Clear the residual of the node, to avoid adding on previous values ---*/
            node[iPoint]->Clear_SurfaceLoad_Res();
          }
        }
        break;
      case DAMPER_BOUNDARY:
        /*--- For all the vertices in the marker iMarker ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          /*--- Retrieve the point ID ---*/
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          /*--- Clear the residual of the node, to avoid adding on previous values ---*/
          node[iPoint]->Clear_SurfaceLoad_Res();
        }
        break;
    }
  }
  
}

void CFEM_ElasticitySolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEM_ElasticitySolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  nPoint = geometry[MESH_0]->GetnPoint();
  
  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/
  
  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_OldSolution();
  }
  
  
}

void CFEM_ElasticitySolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  nPoint = geometry[MESH_0]->GetnPoint();
  
  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/
  
  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_Solution();
  }
  
}

void CFEM_ElasticitySolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double *Kab = NULL;
  unsigned short NelNodes, jNode;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- Compute the components of the jacobian and the stress term ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
          }
        }
        
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
        
      }
      
    }
    
  }
  
  
}

void CFEM_ElasticitySolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  su2double Ks_ab;
  su2double *Kab = NULL;
  su2double *Kk_ab = NULL;
  su2double *Ta = NULL;
  
  su2double *Ta_DE = NULL;
  su2double Ks_ab_DE = 0.0;
  
  unsigned short NelNodes, jNode;
  
  bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
  bool de_effects = config->GetDE_Effects();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;

        /*--- Set current coordinate ---*/
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (de_effects) element_container[DE_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (incompressible) element_container[INC_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);

        /*--- Set reference coordinate ---*/
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
          if (de_effects) element_container[DE_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
          if (incompressible) element_container[INC_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
          if (de_effects) element_container[DE_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
          if (incompressible) element_container[INC_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    if (de_effects) element_container[DE_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    if (incompressible) element_container[INC_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- If incompressible, we compute the Mean Dilatation term first so the volume is already computed ---*/
    if (incompressible) numerics[FEA_TERM]->Compute_MeanDilatation_Term(element_container[INC_TERM][EL_KIND], config);
    
    /*--- Compute the components of the Jacobian and the stress term for the material ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    
    /*--- Compute the electric component of the Jacobian and the stress term ---*/
    if (de_effects) numerics[DE_TERM]->Compute_Tangent_Matrix(element_container[DE_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
      
      /*--- Check if this is my node or not ---*/
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
      /*--- Retrieve the electric contribution to the Residual ---*/
      if (de_effects){
        Ta_DE = element_container[DE_TERM][EL_KIND]->Get_Kt_a(iNode);
        for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta_DE[iVar];
        LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
        
      }
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        /*--- Retrieve the values of the FEA term ---*/
        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        Ks_ab = element_container[FEA_TERM][EL_KIND]->Get_Ks_ab(iNode,jNode);
        if (incompressible) Kk_ab = element_container[INC_TERM][EL_KIND]->Get_Kk_ab(iNode,jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_s_ij[iVar][iVar] = Ks_ab;
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_c_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
            if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
          }
        }
        
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
        if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);
        
        /*--- Retrieve the electric contribution to the Jacobian ---*/
        if (de_effects){
          //          Kab_DE = element_container[DE_TERM][EL_KIND]->Get_Kab(iNode, jNode);
          Ks_ab_DE = element_container[DE_TERM][EL_KIND]->Get_Ks_ab(iNode,jNode);
          
          for (iVar = 0; iVar < nVar; iVar++){
            Jacobian_s_ij[iVar][iVar] = Ks_ab_DE;
          }
          
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
        }
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double Mab;
  unsigned short NelNodes, jNode;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}
    
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    numerics[FEA_TERM]->Compute_Mass_Matrix(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Mab = element_container[FEA_TERM][EL_KIND]->Get_Mab(iNode, jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          MassMatrix_ij[iVar][iVar] = Mab;
        }
        
        MassMatrix.AddBlock(indexNode[iNode], indexNode[jNode], MassMatrix_ij);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Compute_MassRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {

  unsigned long iElem, iVar, iPoint;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;

  su2double Mab;
  unsigned short NelNodes, jNode;

  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    TimeRes.SetBlock_Zero(iPoint);
  }

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}

    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }

    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);

    numerics[FEA_TERM]->Compute_Mass_Matrix(element_container[FEA_TERM][EL_KIND], config);

    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {

      for (jNode = 0; jNode < NelNodes; jNode++) {

        Mab = element_container[FEA_TERM][EL_KIND]->Get_Mab(iNode, jNode);

        for (iVar = 0; iVar < nVar; iVar++) {
          Residual_i[iVar] = Mab * TimeRes_Aux.GetBlock(indexNode[iNode],iVar);
          Residual_j[iVar] = Mab * TimeRes_Aux.GetBlock(indexNode[jNode],iVar);
        }

        TimeRes.AddBlock(indexNode[iNode],Residual_i);
        TimeRes.AddBlock(indexNode[jNode],Residual_j);

      }

    }

  }

}

void CFEM_ElasticitySolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  su2double *Ta = NULL;
  unsigned short NelNodes;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- Compute the components of the jacobian and the stress term ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_NodalStress_Term(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_NodalStress_Term(element_container[FEA_TERM][EL_KIND], config);
    }
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
      
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
    }
    
  }

}

void CFEM_ElasticitySolver::Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iPoint, iElem, iVar;
  unsigned short iNode, iDim, iStress;
  unsigned short nNodes = 0, nStress;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  if (nDim == 2) nStress = 3;
  else nStress = 6;
  
  su2double *Ta = NULL;
  
  unsigned short NelNodes;
  
  /*--- Restart stress to avoid adding results from previous time steps ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iStress = 0; iStress < nStress; iStress++) {
      node[iPoint]->SetStress_FEM(iStress, 0.0);
    }
  }
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      //      for (iDim = 0; iDim < nDim; iDim++) {
      //        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
      //        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
      //        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      //        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      //      }
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- Compute the components of the jacobian and the stress term ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_Averaged_NodalStress(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_Averaged_NodalStress(element_container[FEA_TERM][EL_KIND], config);
    }
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      /*--- This only works if the problem is nonlinear ---*/
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
      
      LinSysReact.AddBlock(indexNode[iNode], Res_Stress_i);
      
      for (iStress = 0; iStress < nStress; iStress++) {
        node[indexNode[iNode]]->AddStress_FEM(iStress,
                                              (element_container[FEA_TERM][EL_KIND]->Get_NodalStress(iNode, iStress) /
                                               geometry->node[indexNode[iNode]]->GetnElem()) );
      }
      
    }
    
  }
  
  su2double *Stress;
  su2double VonMises_Stress, MaxVonMises_Stress = 0.0;
  su2double Sxx,Syy,Szz,Sxy,Sxz,Syz,S1,S2;
  
  /*--- For the number of nodes in the mesh ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Get the stresses, added up from all the elements that connect to the node ---*/
    
    Stress  = node[iPoint]->GetStress_FEM();
    
    /*--- Compute the stress averaged from all the elements connecting to the node and the Von Mises stress ---*/
    
    if (nDim == 2) {
      
      Sxx=Stress[0];
      Syy=Stress[1];
      Sxy=Stress[2];
      
      S1=(Sxx+Syy)/2+sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      S2=(Sxx+Syy)/2-sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      
      VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);
      
    }
    else {
      
      Sxx = Stress[0];
      Syy = Stress[1];
      Szz = Stress[3];
      
      Sxy = Stress[2];
      Sxz = Stress[4];
      Syz = Stress[5];
      
      VonMises_Stress = sqrt(0.5*(   pow(Sxx - Syy, 2.0)
                                  + pow(Syy - Szz, 2.0)
                                  + pow(Szz - Sxx, 2.0)
                                  + 6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)
                                  ));
      
    }
    
    node[iPoint]->SetVonMises_Stress(VonMises_Stress);
    
    /*--- Compute the maximum value of the Von Mises Stress ---*/
    
    MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);
    
  }
  
#ifdef HAVE_MPI
  
  /*--- Compute MaxVonMises_Stress using all the nodes ---*/
  
  su2double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
  SU2_MPI::Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
#endif
  
  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/
  
  Total_CFEA = MaxVonMises_Stress;
  
  
  bool outputReactions = false;
  
  if (outputReactions) {
    
    ofstream myfile;
    myfile.open ("Reactions.txt");
    
    unsigned short iMarker;
    unsigned long iVertex;
    su2double val_Reaction;
    
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
    
    if (!dynamic) {
      /*--- Loop over all the markers  ---*/
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        switch (config->GetMarker_All_KindBC(iMarker)) {
            
            /*--- If it corresponds to a clamped boundary  ---*/
            
          case CLAMPED_BOUNDARY:
            
            myfile << "MARKER " << iMarker << ":" << endl;
            
            /*--- Loop over all the vertices  ---*/
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              
              /*--- Get node index ---*/
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              
              myfile << "Node " << iPoint << "." << " \t ";
              
              for (iDim = 0; iDim < nDim; iDim++) {
                /*--- Retrieve coordinate ---*/
                val_Coord = geometry->node[iPoint]->GetCoord(iDim);
                myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
              }
              
              for (iVar = 0; iVar < nVar; iVar++) {
                /*--- Retrieve reaction ---*/
                val_Reaction = LinSysReact.GetBlock(iPoint, iVar);
                myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
              }
              
              myfile << endl;
            }
            myfile << endl;
            break;
        }
    }
    else if (dynamic) {
      
      switch (config->GetKind_TimeIntScheme_FEA()) {
        case (CD_EXPLICIT):
          cout << "NOT IMPLEMENTED YET" << endl;
          break;
        case (NEWMARK_IMPLICIT):
          
          /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
          if (linear_analysis) {
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
              for (iVar = 0; iVar < nVar; iVar++) {
                Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+    //a0*U(t)
                a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+  //a2*U'(t)
                a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
              }
              TimeRes_Aux.SetBlock(iPoint, Residual);
            }
          }
          else if (nonlinear_analysis) {
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
              for (iVar = 0; iVar < nVar; iVar++) {
                Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
                - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
                + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
                + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
              }
              TimeRes_Aux.SetBlock(iPoint, Residual);
            }
          }
          /*--- Once computed, compute M*TimeRes_Aux ---*/
          MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
          
          /*--- Loop over all the markers  ---*/
          for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            switch (config->GetMarker_All_KindBC(iMarker)) {
                
                /*--- If it corresponds to a clamped boundary  ---*/
                
              case CLAMPED_BOUNDARY:
                
                myfile << "MARKER " << iMarker << ":" << endl;
                
                /*--- Loop over all the vertices  ---*/
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                  
                  /*--- Get node index ---*/
                  iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                  
                  myfile << "Node " << iPoint << "." << " \t ";
                  
                  for (iDim = 0; iDim < nDim; iDim++) {
                    /*--- Retrieve coordinate ---*/
                    val_Coord = geometry->node[iPoint]->GetCoord(iDim);
                    myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
                  }
                  
                  for (iVar = 0; iVar < nVar; iVar++) {
                    /*--- Retrieve the time contribution ---*/
                    Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
                    /*--- Retrieve reaction ---*/
                    val_Reaction = LinSysReact.GetBlock(iPoint, iVar) + Res_Time_Cont[iVar];
                    myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
                  }
                  
                  myfile << endl;
                }
                myfile << endl;
                break;
            }
          
          
          break;
        case (GENERALIZED_ALPHA):
          cout << "NOT IMPLEMENTED YET" << endl;
          break;
      }
      
    }
    
    
    
    myfile.close();
    
  }
  
}

void CFEM_ElasticitySolver::Compute_DeadLoad(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double *Dead_Load = NULL;
  unsigned short NelNodes;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    numerics[FEA_TERM]->Compute_Dead_Load(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Dead_Load = element_container[FEA_TERM][EL_KIND]->Get_FDL_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Dead_Load[iVar] = Dead_Load[iVar];
      
      node[indexNode[iNode]]->Add_BodyForces_Res(Res_Dead_Load);
      
    }
    
  }
  
  
}

void CFEM_ElasticitySolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
}

void CFEM_ElasticitySolver::Compute_IntegrationConstants(CConfig *config) {
  
  su2double Delta_t= config->GetDelta_DynTime();
  
  su2double delta = config->GetNewmark_delta(), alpha = config->GetNewmark_alpha();
  
  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      cout << "NOT IMPLEMENTED YET" << endl;
      break;
    case (NEWMARK_IMPLICIT):
      
      /*--- Integration constants for Newmark scheme ---*/
      
      a_dt[0]= 1 / (alpha*pow(Delta_t,2.0));
      a_dt[1]= delta / (alpha*Delta_t);
      a_dt[2]= 1 / (alpha*Delta_t);
      a_dt[3]= 1 /(2*alpha) - 1;
      a_dt[4]= delta/alpha - 1;
      a_dt[5]= (Delta_t/2) * (delta/alpha - 2);
      a_dt[6]= Delta_t * (1-delta);
      a_dt[7]= delta * Delta_t;
      a_dt[8]= 0.0;
      
      break;
      
    case (GENERALIZED_ALPHA):
      
      /*--- Integration constants for Generalized Alpha ---*/
      /*--- Needs to be updated if accounting for structural damping ---*/
      
      //      su2double beta = config->Get_Int_Coeffs(0);
      //      //  su2double gamma =  config->Get_Int_Coeffs(1);
      //      su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
      //
      //      a_dt[0]= (1 / (beta*pow(Delta_t,2.0))) * ((1 - alpha_m) / (1 - alpha_f)) ;
      //      a_dt[1]= 0.0 ;
      //      a_dt[2]= (1 - alpha_m) / (beta*Delta_t);
      //      a_dt[3]= ((1 - 2*beta)*(1-alpha_m) / (2*beta)) - alpha_m;
      //      a_dt[4]= 0.0;
      //      a_dt[5]= 0.0;
      //      a_dt[6]= Delta_t * (1-delta);
      //      a_dt[7]= delta * Delta_t;
      //      a_dt[8]= (1 - alpha_m) / (beta*pow(Delta_t,2.0));
      
      break;
  }
  
  
}


void CFEM_ElasticitySolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned long iVar, jVar;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node index ---*/
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      if (nDim == 2) {
        Solution[0] = 0.0;  Solution[1] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;
      }
      else {
        Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
      }
      
      node[iPoint]->SetSolution(Solution);
      
      if (dynamic) {
        node[iPoint]->SetSolution_Vel(Solution);
        node[iPoint]->SetSolution_Accel(Solution);
      }
      
      
      /*--- Initialize the reaction vector ---*/
      LinSysReact.SetBlock(iPoint, Residual);
      
      LinSysRes.SetBlock(iPoint, Residual);
      LinSysSol.SetBlock(iPoint, Solution);
      
      /*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/
      
      /*--- Delete the columns for a particular node ---*/
      
      for (iVar = 0; iVar < nPoint; iVar++) {
        if (iVar==iPoint) {
          Jacobian.SetBlock(iVar,iPoint,mId_Aux);
        }
        else {
          Jacobian.SetBlock(iVar,iPoint,mZeros_Aux);
        }
      }
      
      /*--- Delete the rows for a particular node ---*/
      for (jVar = 0; jVar < nPoint; jVar++) {
        if (iPoint!=jVar) {
          Jacobian.SetBlock(iPoint,jVar,mZeros_Aux);
        }
      }
      
      /*--- If the problem is dynamic ---*/
      /*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/
      
      if(dynamic) {
        
        node[iPoint]->SetSolution_time_n(Solution);
        node[iPoint]->SetSolution_Vel_time_n(Solution);
        node[iPoint]->SetSolution_Accel_time_n(Solution);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                            unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node index ---*/
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (nDim == 2) {
      Solution[0] = 0.0;  Solution[1] = 0.0;
    }
    else {
      Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
    }
    
    node[iPoint]->SetSolution(Solution);
    
    if (dynamic) {
      node[iPoint]->SetSolution_Vel(Solution);
      node[iPoint]->SetSolution_Accel(Solution);
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_DispDir(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                            unsigned short val_marker) {

  unsigned short iDim, jDim;

  su2double DispDirVal = config->GetDisp_Dir_Value(config->GetMarker_All_TagBound(val_marker));
  su2double DispDirMult = config->GetDisp_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
  su2double *Disp_Dir_Local= config->GetDisp_Dir(config->GetMarker_All_TagBound(val_marker));
  su2double Disp_Dir_Unit[3] = {0.0, 0.0, 0.0}, Disp_Dir[3] = {0.0, 0.0, 0.0};
  su2double Disp_Dir_Mod = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    Disp_Dir_Mod += Disp_Dir_Local[iDim]*Disp_Dir_Local[iDim];

  Disp_Dir_Mod = sqrt(Disp_Dir_Mod);

  for (iDim = 0; iDim < nDim; iDim++)
    Disp_Dir_Unit[iDim] = Disp_Dir_Local[iDim] / Disp_Dir_Mod;

  su2double TotalDisp;

  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl;

  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  su2double Transfer_Time = 0.0;

  if (Ramp_Load){
    if (Ramp_Time == 0.0)
      ModAmpl = 1.0;
    else
      Transfer_Time = CurrentTime / Ramp_Time;

    switch (config->GetDynamic_LoadTransfer()) {
    case INSTANTANEOUS:
      ModAmpl = 1.0;
      break;
    case POL_ORDER_1:
      ModAmpl = Transfer_Time;
      break;
    case POL_ORDER_3:
      ModAmpl = -2.0 * pow(Transfer_Time,3.0) + 3.0 * pow(Transfer_Time,2.0);
      break;
    case POL_ORDER_5:
      ModAmpl = 6.0 * pow(Transfer_Time, 5.0) - 15.0 * pow(Transfer_Time, 4.0) + 10 * pow(Transfer_Time, 3.0);
      break;
    case SIGMOID_10:
      ModAmpl = (1 / (1+exp(-1.0 * 10.0 * (Transfer_Time - 0.5)) ) );
      break;
    case SIGMOID_20:
      ModAmpl = (1 / (1+exp(-1.0 * 20.0 * (Transfer_Time - 0.5)) ) );
      break;
    }
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else{
    ModAmpl = 1.0;
  }

  TotalDisp = ModAmpl * DispDirVal * DispDirMult;

  for (iDim = 0; iDim < nDim; iDim++)
    Disp_Dir[iDim] = TotalDisp * Disp_Dir_Unit[iDim];

  unsigned long iNode, iVertex;
  unsigned long iPoint, jPoint;

  su2double valJacobian_ij_00 = 0.0;
  su2double auxJacobian_ij[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iNode]->GetDomain()) {

      if (nDim == 2) {
        Solution[0] = Disp_Dir[0];  Solution[1] = Disp_Dir[1];
        Residual[0] = Disp_Dir[0];  Residual[1] = Disp_Dir[1];
      }
      else {
        Solution[0] = Disp_Dir[0];  Solution[1] = Disp_Dir[1];  Solution[2] = Disp_Dir[2];
        Residual[0] = Disp_Dir[0];  Residual[1] = Disp_Dir[1];  Residual[2] = Disp_Dir[2];
      }

      /*--- Initialize the reaction vector ---*/

      LinSysRes.SetBlock(iNode, Residual);
      LinSysSol.SetBlock(iNode, Solution);

      /*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/

      /*--- Delete the full row for node iNode ---*/
      for (jPoint = 0; jPoint < nPoint; jPoint++){

        /*--- Check whether the block is non-zero ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iNode, jPoint,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          if (iNode != jPoint) {
            Jacobian.SetBlock(iNode,jPoint,mZeros_Aux);
          }
          else{
            Jacobian.SetBlock(iNode,jPoint,mId_Aux);
          }
        }
      }

      /*--- Delete the columns for a particular node ---*/

      for (iPoint = 0; iPoint < nPoint; iPoint++){

        /*--- Check if the term K(iPoint, iNode) is 0 ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iPoint,iNode,0,0);

        /*--- If the node iNode has a crossed dependency with the point iPoint ---*/
        if (valJacobian_ij_00 != 0.0 ){

          /*--- Retrieve the Jacobian term ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            for (jDim = 0; jDim < nDim; jDim++){
              auxJacobian_ij[iDim][jDim] = Jacobian.GetBlock(iPoint,iNode,iDim,jDim);
            }
          }

          /*--- Multiply by the imposed displacement ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            Residual[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++){
              Residual[iDim] += auxJacobian_ij[iDim][jDim] * Disp_Dir[jDim];
            }
          }

          if (iNode != iPoint) {
            /*--- The term is substracted from the residual (right hand side) ---*/
            LinSysRes.SubtractBlock(iPoint, Residual);
            /*--- The Jacobian term is now set to 0 ---*/
            Jacobian.SetBlock(iPoint,iNode,mZeros_Aux);
          }

        }

      }

    }

  }


}

void CFEM_ElasticitySolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                                           unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  bool first_iter = (config->GetIntIter() == 0);
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);    // Nonlinear analysis.
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  
  su2double solNorm = 0.0, solNorm_recv = 0.0;
  
#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (disc_adj_fem) {

    if (nonlinear_analysis) {

      /*--- For nonlinear discrete adjoint, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)): ABSOLUTE, norm of the incremental displacements ------------*/
      /*--- RTOL = norm(Residual(k): ABSOLUTE, norm of the residual (T-F) -----------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k): ABSOLUTE, norm of the product disp * res -----------*/

      Conv_Check[0] = LinSysSol.norm();               // Norm of the delta-solution vector
      Conv_Check[1] = LinSysRes.norm();               // Norm of the residual
      Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);
    }
    else {
      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/
      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);

    }

  }
  else {
    if (nonlinear_analysis){

      /*--- If the problem is nonlinear, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
      /*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/

      if (first_iter){
        Conv_Ref[0] = 1.0;                                        // Position for the norm of the solution
        Conv_Ref[1] = max(LinSysRes.norm(), EPS);                 // Position for the norm of the residual
        Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);    // Position for the energy tolerance

        /*--- Make sure the computation runs at least 2 iterations ---*/
        Conv_Check[0] = 1.0;
        Conv_Check[1] = 1.0;
        Conv_Check[2] = 1.0;
      }
      else {
        /*--- Compute the norm of the solution vector Uk ---*/
        for (iPoint = 0; iPoint < nPointDomain; iPoint++){
          for (iVar = 0; iVar < nVar; iVar++){
            solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
          }
        }

        // We need to communicate the norm of the solution and compute the RMS throughout the different processors

#ifdef HAVE_MPI
        /*--- We sum the squares of the norms across the different processors ---*/
        SU2_MPI::Allreduce(&solNorm, &solNorm_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        solNorm_recv         = solNorm;
#endif

        Conv_Ref[0] = max(sqrt(solNorm_recv), EPS);           // Norm of the solution vector

        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_RELATIVE:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
          default:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
        }

      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

    } else {

      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                                   unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                           unsigned short val_marker) {
  
  /*--- Retrieve the normal pressure and the application conditions for the considered boundary ---*/
  
  su2double NormalLoad = config->GetLoad_Value(config->GetMarker_All_TagBound(val_marker));
  su2double TotalLoad = 0.0;
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl;
  
  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  su2double Transfer_Time = 0.0;
  
  if (Ramp_Load) {
    if (Ramp_Time == 0.0)
      ModAmpl = 1.0;
    else
      Transfer_Time = CurrentTime / Ramp_Time;
    
    switch (config->GetDynamic_LoadTransfer()) {
    case INSTANTANEOUS:
      ModAmpl = 1.0;
      break;
    case POL_ORDER_1:
      ModAmpl = Transfer_Time;
      break;
    case POL_ORDER_3:
      ModAmpl = -2.0 * pow(Transfer_Time,3.0) + 3.0 * pow(Transfer_Time,2.0);
      if (CurrentTime > Ramp_Time) ModAmpl = 1.0;
      break;
    case POL_ORDER_5:
      ModAmpl = 6.0 * pow(Transfer_Time, 5.0) - 15.0 * pow(Transfer_Time, 4.0) + 10 * pow(Transfer_Time, 3.0);
      if (CurrentTime > Ramp_Time) ModAmpl = 1.0;
      break;
    case SIGMOID_10:
      ModAmpl = (1 / (1+exp(-1.0 * 10.0 * (Transfer_Time - 0.5)) ) );
      break;
    case SIGMOID_20:
      ModAmpl = (1 / (1+exp(-1.0 * 20.0 * (Transfer_Time - 0.5)) ) );
      break;
    }
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else {
    ModAmpl = 1.0;
  }

  TotalLoad = ModAmpl * NormalLoad;
  
  /*--- Do only if there is a load applied.
   *--- This reduces the computational cost for cases in which we want boundaries with no load.
   */
  
  if (TotalLoad != 0.0) {
    
    unsigned long iElem;
    unsigned short nNodes = 0;
    su2double Length_Elem_ref = 0.0,  Area_Elem_ref = 0.0;
    su2double Length_Elem_curr = 0.0, Area_Elem_curr = 0.0;
    unsigned long indexNode[4]   = {0,0,0,0};
    unsigned long indexVertex[4] = {0,0,0,0};
    su2double nodeCoord_ref[4][3], nodeCoord_curr[4][3];
    su2double *nodal_normal, nodal_normal_unit[3];
    su2double normal_ref_unit[3], normal_curr_unit[3];
    su2double Norm, dot_Prod;
    su2double val_Coord, val_Sol;
    unsigned short iNode, iDim;
    unsigned long iVertex, iPoint;
    su2double a_ref[3], b_ref[3], AC_ref[3], BD_ref[3];
    su2double a_curr[3], b_curr[3], AC_curr[3], BD_curr[3];
    
    /*--- Determine whether the load conditions are applied in the reference or in the current configuration ---*/
    
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS); // Nonlinear analysis.
    
    for (iNode = 0; iNode < 4; iNode++) {
      for (iDim = 0; iDim < 3; iDim++) {
        nodeCoord_ref[iNode][iDim]  = 0.0;
        nodeCoord_curr[iNode][iDim] = 0.0;
      }
    }
    
    for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
      
      /*--- Identify the kind of boundary element ---*/
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == LINE)           nNodes = 2;
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       nNodes = 3;
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL)  nNodes = 4;
      
      /*--- Retrieve the boundary reference and current coordinates ---*/
      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; iDim++) {
          val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
          val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
          /*--- Assign values to the container ---*/
          nodeCoord_ref[iNode][iDim]  = val_Coord;
          nodeCoord_curr[iNode][iDim] = val_Sol;
        }
      }
      
      /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        for (iNode = 0; iNode < nNodes; iNode++) {
          if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
        }
      }
      
      /*--- Retrieve the reference normal for one of the points. They go INSIDE the structural domain. ---*/
      nodal_normal = geometry->vertex[val_marker][indexVertex[0]]->GetNormal();
      Norm = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Norm += nodal_normal[iDim]*nodal_normal[iDim];
      }
      Norm = sqrt(Norm);
      for (iDim = 0; iDim < nDim; iDim++) {
        nodal_normal_unit[iDim] = nodal_normal[iDim] / Norm;
      }
      
      /*--- Compute area (3D), and length of the surfaces (2D), and the unitary normal vector in current configuration ---*/
      
      if (nDim == 2) {
        
        /*-- Compute the vector a in reference and current configurations ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          a_ref[iDim]  = nodeCoord_ref[0][iDim] -nodeCoord_ref[1][iDim];
          a_curr[iDim] = nodeCoord_curr[0][iDim]-nodeCoord_curr[1][iDim];
        }
        
        /*-- Compute the length of the boundary element in reference and current configurations ---*/
        Length_Elem_curr = sqrt(a_curr[0]*a_curr[0]+a_curr[1]*a_curr[1]);
        Length_Elem_ref  = sqrt(a_ref[0]*a_ref[0]+a_ref[1]*a_ref[1]);
        
        /*-- Compute the length of the boundary element in reference and current configurations ---*/
        normal_ref_unit[0] =   a_ref[1] /Length_Elem_ref;
        normal_ref_unit[1] = -(a_ref[0])/Length_Elem_ref;
        
        normal_curr_unit[0] =   a_curr[1] /Length_Elem_curr;
        normal_curr_unit[1] = -(a_curr[0])/Length_Elem_curr;
        
        /*-- Dot product to check the element orientation in the reference configuration ---*/
        dot_Prod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
        }
        
        /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
        /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
        /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
        if (dot_Prod < 0) {
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
            normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
          }
        }
        
        if (linear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_ref * normal_ref_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_ref * normal_ref_unit[1];
          
          node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
          node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
        }
        else if (nonlinear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[1];
          
          node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
          node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
        }
        
      }
      
      if (nDim == 3) {
        
        if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
          
          for (iDim = 0; iDim < nDim; iDim++) {
            a_ref[iDim] = nodeCoord_ref[1][iDim]-nodeCoord_ref[0][iDim];
            b_ref[iDim] = nodeCoord_ref[2][iDim]-nodeCoord_ref[0][iDim];
            
            a_curr[iDim] = nodeCoord_curr[1][iDim]-nodeCoord_curr[0][iDim];
            b_curr[iDim] = nodeCoord_curr[2][iDim]-nodeCoord_curr[0][iDim];
          }
          
          su2double Ni=0, Nj=0, Nk=0;
          
          /*--- Reference configuration ---*/
          Ni = a_ref[1]*b_ref[2] - a_ref[2]*b_ref[1];
          Nj = a_ref[2]*b_ref[0] - a_ref[0]*b_ref[2];
          Nk = a_ref[0]*b_ref[1] - a_ref[1]*b_ref[0];
          
          Area_Elem_ref = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_ref_unit[0] = Ni / Area_Elem_ref;
          normal_ref_unit[1] = Nj / Area_Elem_ref;
          normal_ref_unit[2] = Nk / Area_Elem_ref;
          
          /*--- Current configuration ---*/
          Ni = a_curr[1]*b_curr[2] - a_curr[2]*b_curr[1];
          Nj = a_curr[2]*b_curr[0] - a_curr[0]*b_curr[2];
          Nk = a_curr[0]*b_curr[1] - a_curr[1]*b_curr[0];
          
          Area_Elem_curr = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_curr_unit[0] = Ni / Area_Elem_curr;
          normal_curr_unit[1] = Nj / Area_Elem_curr;
          normal_curr_unit[2] = Nk / Area_Elem_curr;
          
          /*-- Dot product to check the element orientation in the reference configuration ---*/
          dot_Prod = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
          }
          
          /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
          /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
          /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
          if (dot_Prod < 0) {
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
              normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
            }
          }
          
          if (linear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
          }
          
        }
        
        else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
          
          for (iDim = 0; iDim < nDim; iDim++) {
            AC_ref[iDim] = nodeCoord_ref[2][iDim]-nodeCoord_ref[0][iDim];
            BD_ref[iDim] = nodeCoord_ref[3][iDim]-nodeCoord_ref[1][iDim];
            
            AC_curr[iDim] = nodeCoord_curr[2][iDim]-nodeCoord_curr[0][iDim];
            BD_curr[iDim] = nodeCoord_curr[3][iDim]-nodeCoord_curr[1][iDim];
          }
          
          su2double Ni=0, Nj=0, Nk=0;
          
          /*--- Reference configuration ---*/
          Ni=AC_ref[1]*BD_ref[2]-AC_ref[2]*BD_ref[1];
          Nj=-AC_ref[0]*BD_ref[2]+AC_ref[2]*BD_ref[0];
          Nk=AC_ref[0]*BD_ref[1]-AC_ref[1]*BD_ref[0];
          
          Area_Elem_ref = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_ref_unit[0] = Ni / Area_Elem_ref;
          normal_ref_unit[1] = Nj / Area_Elem_ref;
          normal_ref_unit[2] = Nk / Area_Elem_ref;
          
          /*--- Current configuration ---*/
          Ni=AC_curr[1]*BD_curr[2]-AC_curr[2]*BD_curr[1];
          Nj=-AC_curr[0]*BD_curr[2]+AC_curr[2]*BD_curr[0];
          Nk=AC_curr[0]*BD_curr[1]-AC_curr[1]*BD_curr[0];
          
          Area_Elem_curr = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_curr_unit[0] = Ni / Area_Elem_curr;
          normal_curr_unit[1] = Nj / Area_Elem_curr;
          normal_curr_unit[2] = Nk / Area_Elem_curr;
          
          /*-- Dot product to check the element orientation in the reference configuration ---*/
          dot_Prod = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
          }
          
          /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
          /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
          /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
          if (dot_Prod < 0) {
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
              normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
            }
          }
          
          if (linear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[3]]->Add_SurfaceLoad_Res(Residual);
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[3]]->Add_SurfaceLoad_Res(Residual);
          }
          
        }
        
      }
      
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  su2double a[3], b[3], AC[3], BD[3];
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
  su2double Length_Elem = 0.0, Area_Elem = 0.0;
  unsigned short iDim;
  
  su2double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
  su2double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
  su2double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));
  
  su2double TotalLoad;
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl;
  
  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  su2double Transfer_Time = 0.0;

  if (Ramp_Load) {
    if (Ramp_Time == 0.0)
      ModAmpl = 1.0;
    else
      Transfer_Time = CurrentTime / Ramp_Time;

    switch (config->GetDynamic_LoadTransfer()) {
    case INSTANTANEOUS:
      ModAmpl = 1.0;
      break;
    case POL_ORDER_1:
      ModAmpl = Transfer_Time;
      break;
    case POL_ORDER_3:
      ModAmpl = -2.0 * pow(Transfer_Time,3.0) + 3.0 * pow(Transfer_Time,2.0);
      if (Transfer_Time > 1.0) ModAmpl = 1.0;
      break;
    case POL_ORDER_5:
      ModAmpl = 6.0 * pow(Transfer_Time, 5.0) - 15.0 * pow(Transfer_Time, 4.0) + 10 * pow(Transfer_Time, 3.0);
      if (Transfer_Time > 1.0) ModAmpl = 1.0;
      break;
    case SIGMOID_10:
      ModAmpl = (1 / (1+exp(-1.0 * 10.0 * (Transfer_Time - 0.5)) ) );
      break;
    case SIGMOID_20:
      ModAmpl = (1 / (1+exp(-1.0 * 20.0 * (Transfer_Time - 0.5)) ) );
      break;
    }
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else {
    ModAmpl = 1.0;
  }

  TotalLoad = ModAmpl * LoadDirVal * LoadDirMult;
  
  /*--- Compute the norm of the vector that was passed in the config file ---*/
  su2double Norm = 1.0;
  if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
  if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);
  
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
    if (nDim == 3) {
      
      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);  Coord_3 = geometry->node[Point_3]->GetCoord();
      }
      
    }
    
    /*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
    if (nDim == 2) {
      
      for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];
      
      Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
      //      Normal_Elem[0] =   a[1];
      //      Normal_Elem[1] = -(a[0]);
      
    }
    
    if (nDim == 3) {
      
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
        
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_1[iDim]-Coord_0[iDim];
          b[iDim] = Coord_2[iDim]-Coord_0[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=a[1]*b[2]-a[2]*b[1];
        Nj=-a[0]*b[2]+a[2]*b[0];
        Nk=a[0]*b[1]-a[1]*b[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
      
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        
        for (iDim = 0; iDim < nDim; iDim++) {
          AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
          BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=AC[1]*BD[2]-AC[2]*BD[1];
        Nj=-AC[0]*BD[2]+AC[2]*BD[0];
        Nk=AC[0]*BD[1]-AC[1]*BD[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
    }
    
    if (nDim == 2) {
      
      Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
      Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
      
      node[Point_0]->Add_SurfaceLoad_Res(Residual);
      node[Point_1]->Add_SurfaceLoad_Res(Residual);
      
    }
    
    else {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
        
        Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);
        node[Point_2]->Add_SurfaceLoad_Res(Residual);
        
      }
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        
        Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);
        node[Point_2]->Add_SurfaceLoad_Res(Residual);
        node[Point_3]->Add_SurfaceLoad_Res(Residual);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                         unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Damper(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                      unsigned short val_marker) {

  unsigned short iVar;
  su2double dampValue, dampC;
  su2double dampConstant = config->GetDamper_Constant(config->GetMarker_All_TagBound(val_marker));
  unsigned long nVertex = geometry->GetnVertex(val_marker);

  /*--- The damping is distributed evenly over all the nodes in the marker ---*/
  dampC = dampConstant / (nVertex + EPS);

  unsigned long Point_0, Point_1, Point_2, Point_3;
  unsigned long iElem;

  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);

    for (iVar = 0; iVar < nVar; iVar++){

        dampValue = - 1.0 * dampC * node[Point_0]->GetSolution_Vel(iVar);
        node[Point_0]->Set_SurfaceLoad_Res(iVar, dampValue);

        dampValue = - 1.0 * dampC * node[Point_1]->GetSolution_Vel(iVar);
        node[Point_1]->Set_SurfaceLoad_Res(iVar, dampValue);
    }

    if (nDim == 3) {

      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);

      for (iVar = 0; iVar < nVar; iVar++){
          dampValue = - 1.0 * dampC * node[Point_2]->GetSolution_Vel(iVar);
          node[Point_2]->Set_SurfaceLoad_Res(iVar, dampValue);
      }

      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);
        for (iVar = 0; iVar < nVar; iVar++){
            dampValue = - 1.0 * dampC * node[Point_3]->GetSolution_Vel(iVar);
            node[Point_3]->Set_SurfaceLoad_Res(iVar, dampValue);
        }
      }

    }

  }


}

void CFEM_ElasticitySolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) { }

void CFEM_ElasticitySolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CFEM_ElasticitySolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, jPoint;
  unsigned short iVar, jVar;
  
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation.
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  bool restart = config->GetRestart();                          // Restart solution
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter());  // Restart iteration
  
  bool incremental_load = config->GetIncrementalLoad();
  
  if (!dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Add the external contribution to the residual    ---*/
      /*--- (the terms that are constant over the time step) ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }

      /*---  Add the contribution to the residual due to flow loads (FSI contribution) ---*/
      if (fsi) {
        if (incremental_load){
          for (iVar = 0; iVar < nVar; iVar++){
            Res_FSI_Cont[iVar] = loadIncrement * node[iPoint]->Get_FlowTraction(iVar);
          }
        }
        else {
            for (iVar = 0; iVar < nVar; iVar++){
              Res_FSI_Cont[iVar] = node[iPoint]->Get_FlowTraction(iVar);
            }
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
    
  }
  
  if (dynamic) {
    
    /*--- Add the mass matrix contribution to the Jacobian ---*/
    
    /*
     * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
     * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
     * iteration, as the Jacobian is recomputed
     *
     * If the problem is linear, we add the Mass Matrix contribution to the Jacobian at the first calculation.
     * From then on, the Jacobian is always the same matrix.
     *
     */
    
    if ((nonlinear_analysis && (newton_raphson || first_iter)) ||
        (linear_analysis && initial_calc) ||
        (linear_analysis && restart && initial_calc_restart)) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (jPoint = 0; jPoint < nPoint; jPoint++) {
          for(iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
            }
          }
          Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
        }
      }
    }
    
    
    /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
    if (linear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+    //a0*U(t)
          a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+  //a2*U'(t)
          a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
    }
    else if (nonlinear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
          - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
          + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
          + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
      
    }
    
    /*--- Once computed, compute M*TimeRes_Aux ---*/
    MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
    /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Dynamic contribution ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
      }
      //Res_Time_Cont = TimeRes.GetBlock(iPoint);
      LinSysRes.AddBlock(iPoint, Res_Time_Cont);
      
      /*--- External surface load contribution ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      
      /*--- Body forces contribution (dead load) ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
      /*--- FSI contribution (flow loads) ---*/
      if (fsi) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = loadIncrement * node[iPoint]->Get_FlowTraction(iVar);
          }
        }
        else {
          Res_FSI_Cont = node[iPoint]->Get_FlowTraction();
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
  }
  
  
}

void CFEM_ElasticitySolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);    // Geometrically linear problems
  bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Geometrically non-linear problems
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);          // Dynamic simulations.
  
  /*--- Update solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Displacements component of the solution ---*/
      
      /*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/
      
      if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
      if (nonlinear)  node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
    }
    
  }
  
  if (dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Acceleration component of the solution ---*/
        /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
        
        Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
                                node[iPoint]->GetSolution_time_n(iVar)) -
        a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
        a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      }
      
      /*--- Set the acceleration in the node structure ---*/
      
      node[iPoint]->SetSolution_Accel(Solution);
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Velocity component of the solution ---*/
        /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
        
        Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
        a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
        a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
        
      }
      
      /*--- Set the velocity in the node structure ---*/
      
      node[iPoint]->SetSolution_Vel(Solution);
      
    }
    
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  
}

void CFEM_ElasticitySolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  su2double *valSolutionPred;
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  /*--- Update solution and set it to be the solution after applying relaxation---*/
  
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    
    valSolutionPred = node[iPoint]->GetSolution_Pred();
    
    node[iPoint]->SetSolution(valSolutionPred);
  }
  
  if (dynamic){
    
    /*--- Compute velocities and accelerations ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Acceleration component of the solution ---*/
        /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

        Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
            node[iPoint]->GetSolution_time_n(iVar)) -
            a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
            a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      }

      /*--- Set the acceleration in the node structure ---*/

      node[iPoint]->SetSolution_Accel(Solution);

      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Velocity component of the solution ---*/
        /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

        Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
            a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
            a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);

      }

      /*--- Set the velocity in the node structure ---*/

      node[iPoint]->SetSolution_Vel(Solution);

    }
  
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- After the solution has been communicated, set the 'old' predicted solution as the solution ---*/
  /*--- Loop over n points (as we have already communicated everything ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->SetSolution_Pred_Old(iVar,node[iPoint]->GetSolution(iVar));
    }
  }
  
  
}


void CFEM_ElasticitySolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, jPoint;
  unsigned short iVar, jVar;
  
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation.
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  bool restart = config->GetRestart();                          // Restart solution
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter());  // Restart iteration
  
  su2double alpha_f = config->Get_Int_Coeffs(2);
  
  bool incremental_load = config->GetIncrementalLoad();
  
  if (!dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Add the external contribution to the residual    ---*/
      /*--- (the terms that are constant over the time step) ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
    }
    
  }
  
  if (dynamic) {
    
    /*--- Add the mass matrix contribution to the Jacobian ---*/
    
    /*
     * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
     * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
     * iteration, as the Jacobian is recomputed
     *
     * If the problem is linear, we add the Mass Matrix contribution to the Jacobian at the first calculation.
     * From then on, the Jacobian is always the same matrix.
     *
     */
    
    if ((nonlinear_analysis && (newton_raphson || first_iter)) ||
        (linear_analysis && initial_calc) ||
        (linear_analysis && restart && initial_calc_restart)) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (jPoint = 0; jPoint < nPoint; jPoint++) {
          for(iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
            }
          }
          Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
        }
      }
    }
    
    
    /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
    if (linear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+    //a0*U(t)
          a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+  //a2*U'(t)
          a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
    }
    else if (nonlinear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
          - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
          + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
          + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
    }
    /*--- Once computed, compute M*TimeRes_Aux ---*/
    MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
    /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      /*--- Dynamic contribution ---*/
      //Res_Time_Cont = TimeRes.GetBlock(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
      }
      LinSysRes.AddBlock(iPoint, Res_Time_Cont);
      /*--- External surface load contribution ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * ( (1 - alpha_f) * node[iPoint]->Get_SurfaceLoad_Res(iVar) +
                                                alpha_f  * node[iPoint]->Get_SurfaceLoad_Res_n(iVar) );
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = (1 - alpha_f) * node[iPoint]->Get_SurfaceLoad_Res(iVar) +
          alpha_f  * node[iPoint]->Get_SurfaceLoad_Res_n(iVar);
        }
      }
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces.
       *--- It is constant over time, so it's not necessary to distribute it. ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
      /*--- Add FSI contribution ---*/
      if (fsi) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = loadIncrement * ( (1 - alpha_f) * node[iPoint]->Get_FlowTraction(iVar) +
                                                  alpha_f  * node[iPoint]->Get_FlowTraction_n(iVar) );
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = (1 - alpha_f) * node[iPoint]->Get_FlowTraction(iVar) +
            alpha_f  * node[iPoint]->Get_FlowTraction_n(iVar);
          }
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
  }
  
}

void CFEM_ElasticitySolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);    // Geometrically linear problems
  bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Geometrically non-linear problems
  
  /*--- Update solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Displacements component of the solution ---*/
      
      /*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/
      
      if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
      if (nonlinear)  node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
    }
    
  }
  
  /*--- Perform the MPI communication of the solution, displacements only ---*/
  
  Set_MPI_Solution_DispOnly(geometry, config);
  
}

void CFEM_ElasticitySolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
  
  /*--- Compute solution at t_n+1, and update velocities and accelerations ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U(t+dt) = 1/alpha_f*(U(t+1-alpha_f)-alpha_f*U(t)) ---*/
      
      Solution[iVar]=(1 / (1 - alpha_f))*(node[iPoint]->GetSolution(iVar) -
                                          alpha_f * node[iPoint]->GetSolution_time_n(iVar));
      

    }
    
    /*--- Set the solution in the node structure ---*/
    
    node[iPoint]->SetSolution(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Acceleration component of the solution ---*/
      /*--- U''(t+dt-alpha_m) = a8*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
      
      Solution_Interm[iVar]=a_dt[8]*( node[iPoint]->GetSolution(iVar) -
                                     node[iPoint]->GetSolution_time_n(iVar)) -
      a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
      a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      
      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U''(t+dt) = 1/alpha_m*(U''(t+1-alpha_m)-alpha_m*U''(t)) ---*/
      
      Solution[iVar]=(1 / (1 - alpha_m))*(Solution_Interm[iVar] - alpha_m * node[iPoint]->GetSolution_Accel_time_n(iVar));
    }
    
    /*--- Set the acceleration in the node structure ---*/
    
    node[iPoint]->SetSolution_Accel(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Velocity component of the solution ---*/
      /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
      
      Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
      a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
      a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
      
    }
    
    /*--- Set the velocity in the node structure ---*/
    
    node[iPoint]->SetSolution_Vel(Solution);
    
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
}

void CFEM_ElasticitySolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint;
  bool fsi = config->GetFSI_Simulation();
  
  /*--- Set the load conditions of the time step n+1 as the load conditions for time step n ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->Set_SurfaceLoad_Res_n();
    if (fsi) node[iPoint]->Set_FlowTraction_n();
  }
  
}

void CFEM_ElasticitySolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long IterLinSol = 0, iPoint, total_index;
  unsigned short iVar;
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
    
  }
  
  CSysSolve femSystem;
  IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
}



void CFEM_ElasticitySolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry,
                                        CGeometry **flow_geometry, CConfig *fea_config,
                                        CConfig *flow_config, CNumerics *fea_numerics) {
  
  unsigned short nMarkerFSI, nMarkerFlow;    // Number of markers on FSI problem, FEA and Flow side
  unsigned short iMarkerFSI, iMarkerFlow;    // Variables for iteration over markers
  int Marker_Flow = -1;
  
  unsigned long iVertex, iPoint;                // Variables for iteration over vertices and nodes
  
  unsigned short iDim, jDim;
  
  // Check the kind of fluid problem
  bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
                             (flow_config->GetKind_Solver() == RANS) );
  
  /*--- Redimensionalize the pressure ---*/
  
  su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
  su2double factorForces;
  
  Velocity_Real = flow_config->GetVelocity_FreeStream();
  Density_Real = flow_config->GetDensity_FreeStream();
  
  Velocity_ND = flow_config->GetVelocity_FreeStreamND();
  Density_ND = flow_config->GetDensity_FreeStreamND();
  
  Velocity2_Real = 0.0;
  Velocity2_ND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
    Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
  }
  
  factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);
  
  /*--- Apply a ramp to the transfer of the fluid loads ---*/
  
  su2double ModAmpl;
  su2double CurrentTime = fea_config->GetCurrent_DynTime();
  su2double Static_Time = fea_config->GetStatic_Time();
  
  bool Ramp_Load = fea_config->GetRamp_Load();
  su2double Ramp_Time = fea_config->GetRamp_Time();
  
  if (CurrentTime <= Static_Time) { ModAmpl=0.0; }
  else if((CurrentTime > Static_Time) &&
          (CurrentTime <= (Static_Time + Ramp_Time)) &&
          (Ramp_Load)) {
    ModAmpl = (CurrentTime-Static_Time) / Ramp_Time;
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else { ModAmpl = 1.0; }
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerFSI = (fea_config->GetMarker_n_ZoneInterface())/2;
  nMarkerFlow = flow_geometry[MESH_0]->GetnMarker();    // Retrieve total number of markers on Fluid side
  
  // Parameters for the calculations
  // Pn: Pressure
  // Pinf: Pressure_infinite
  // div_vel: Velocity divergence
  // Dij: Dirac delta
  su2double Pn = 0.0, Pinf = 0.0, div_vel = 0.0, Dij = 0.0;
  su2double Viscosity = 0.0;
  su2double **Grad_PrimVar = NULL;
  su2double Tau[3][3];
  
  unsigned long Point_Flow, Point_Struct;
  su2double *Normal_Flow;
  
  su2double *tn_f;
  tn_f         = new su2double [nVar];      // Fluid traction
  
#ifndef HAVE_MPI
  
  unsigned long nVertexFlow;            // Number of vertices on FEA and Flow side
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->Clear_FlowTraction();
  }
  
  /*--- Loop over all the markers on the interface ---*/
  
  for (iMarkerFSI = 0; iMarkerFSI < nMarkerFSI; iMarkerFSI++) {
    
    /*--- Identification of the markers ---*/
    
    /*--- Current fluid marker ---*/
    for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++) {
      if (flow_config->GetMarker_All_ZoneInterface(iMarkerFlow) == (iMarkerFSI+1)) {
        Marker_Flow = iMarkerFlow;
      }
    }
    
    nVertexFlow = flow_geometry[MESH_0]->GetnVertex(Marker_Flow);  // Retrieve total number of vertices on Fluid marker
    
    /*--- Loop over the nodes in the fluid mesh, calculate the tf vector (unitary) ---*/
    /*--- Here, we are looping over the fluid, and we find the pointer to the structure (Point_Struct) ---*/
    for (iVertex = 0; iVertex < nVertexFlow; iVertex++) {
      
      // Node from the flow mesh
      Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();
      
      // Normals at the vertex: these normals go inside the fluid domain.
      Normal_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNormal();
      
      // Corresponding node on the structural mesh
      Point_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorPoint();
      
      // Retrieve the values of pressure, viscosity and density
      if (incompressible) {
        
        Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
        Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
        
        if (viscous_flow) {
          
          Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
          Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
        }
      }
      else if (compressible) {
        
        Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
        Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
        
        if (viscous_flow) {
          
          Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
          Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
        }
      }
      
      // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
      for (iDim = 0; iDim < nDim; iDim++) {
        tn_f[iDim] = -(Pn-Pinf)*Normal_Flow[iDim];
      }
      
      // Calculate tn in the fluid nodes for the viscous term
      
      if (viscous_flow) {
        
        // Divergence of the velocity
        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
        if (incompressible) div_vel = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          
          for (jDim = 0 ; jDim < nDim; jDim++) {
            // Dirac delta
            Dij = 0.0; if (iDim == jDim) Dij = 1.0;
            
            // Viscous stress
            Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
            TWO3*Viscosity*div_vel*Dij;
            
            // Viscous component in the tn vector --> Units of force (non-dimensional).
            tn_f[iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
          }
        }
      }
      
      // Rescale tn to SI units and apply time-dependent coefficient (static structure, ramp load, full load)
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Residual[iDim] = tn_f[iDim]*factorForces*ModAmpl;
      }
      
      /*--- Set the Flow traction ---*/
      //node[Point_Struct]->Set_FlowTraction(Residual);
      /*--- Add to the Flow traction (to add values to corners...) ---*/
      node[Point_Struct]->Add_FlowTraction(Residual);
    }
    
  }
  
#else
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  unsigned long nLocalVertexStruct = 0, nLocalVertexFlow = 0;
  
  unsigned long MaxLocalVertexStruct = 0, MaxLocalVertexFlow = 0;
  
  unsigned long nBuffer_FlowTraction = 0, nBuffer_StructTraction = 0;
  unsigned long nBuffer_DonorIndices = 0, nBuffer_SetIndex = 0;
  
  unsigned long Processor_Struct;
  
  int iProcessor, nProcessor = 0;
  
  unsigned short nMarkerStruct, iMarkerStruct;    // Variables for iteration over markers
  int Marker_Struct = -1;
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerFSI     = (flow_config->GetMarker_n_ZoneInterface())/2;
  nMarkerStruct  = fea_geometry[MESH_0]->GetnMarker();
  nMarkerFlow    = flow_geometry[MESH_0]->GetnMarker();
  
  nProcessor = size;
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->Clear_FlowTraction();
  }
  
  /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/
  
  for (iMarkerFSI = 1; iMarkerFSI <= nMarkerFSI; iMarkerFSI++) {
    
    Marker_Struct = -1;
    Marker_Flow = -1;
    
    /*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
    unsigned long Buffer_Send_nVertexStruct[1], *Buffer_Recv_nVertexStruct = NULL;
    unsigned long Buffer_Send_nVertexFlow[1], *Buffer_Recv_nVertexFlow = NULL;
    
    /*--- The markers on the fluid and structural side are tagged with the same index.
     *--- This is independent of the MPI domain decomposition.
     *--- We need to loop over all markers on both sides and get the number of nodes
     *--- that belong to each FSI marker for each processor ---*/
    
    /*--- On the structural side ---*/
    
    for (iMarkerStruct = 0; iMarkerStruct < nMarkerStruct; iMarkerStruct++) {
      /*--- If the tag GetMarker_All_ZoneInterface(iMarkerStruct) equals the index we are looping at ---*/
      if ( fea_config->GetMarker_All_ZoneInterface(iMarkerStruct) == iMarkerFSI ) {
        /*--- We have identified the local index of the FEA marker ---*/
        /*--- Store the number of local points that belong to Marker_Struct on each processor ---*/
        /*--- This includes the halo nodes ---*/
        nLocalVertexStruct = fea_geometry[MESH_0]->GetnVertex(iMarkerStruct);
        /*--- Store the identifier for the structural marker ---*/
        Marker_Struct = iMarkerStruct;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
      else {
        /*--- If the tag hasn't matched any tag within the FEA markers ---*/
        nLocalVertexStruct = 0;
        Marker_Struct = -1;
      }
    }
    
    /*--- On the fluid side ---*/
    
    for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++) {
      /*--- If the tag GetMarker_All_ZoneInterface(iMarkerFlow) equals the index we are looping at ---*/
      if ( flow_config->GetMarker_All_ZoneInterface(iMarkerFlow) == iMarkerFSI ) {
        /*--- We have identified the local index of the Flow marker ---*/
        /*--- Store the number of local points that belong to Marker_Flow on each processor ---*/
        /*--- This includes the halo nodes ---*/
        nLocalVertexFlow = flow_geometry[MESH_0]->GetnVertex(iMarkerFlow);
        /*--- Store the identifier for the fluid marker ---*/
        Marker_Flow = iMarkerFlow;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
      else {
        /*--- If the tag hasn't matched any tag within the Flow markers ---*/
        nLocalVertexFlow = 0;
        Marker_Flow = -1;
      }
    }
    
    Buffer_Send_nVertexStruct[0] = nLocalVertexStruct;                  // Retrieve total number of vertices on FEA marker
    Buffer_Send_nVertexFlow[0] = nLocalVertexFlow;                    // Retrieve total number of vertices on Flow marker
    if (rank == MASTER_NODE) Buffer_Recv_nVertexStruct = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
    if (rank == MASTER_NODE) Buffer_Recv_nVertexFlow = new unsigned long[size];     // Allocate memory to receive how many vertices are on each rank on the fluid side
    
    /*--- We receive MaxLocalVertexFEA as the maximum number of vertices in one single processor on the structural side---*/
    SU2_MPI::Allreduce(&nLocalVertexStruct, &MaxLocalVertexStruct, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    /*--- We receive MaxLocalVertexFlow as the maximum number of vertices in one single processor on the fluid side ---*/
    SU2_MPI::Allreduce(&nLocalVertexFlow, &MaxLocalVertexFlow, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    
    /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
    SU2_MPI::Gather(&Buffer_Send_nVertexStruct, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexStruct, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the fluid side ---*/
    SU2_MPI::Gather(&Buffer_Send_nVertexFlow, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexFlow, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- We will be gathering the structural coordinates into the master node ---*/
    /*--- Then we will distribute them using a scatter operation into the appropriate fluid processor ---*/
    nBuffer_FlowTraction = MaxLocalVertexFlow * nDim;
    nBuffer_StructTraction = MaxLocalVertexStruct * nDim;
    
    /*--- We will be gathering donor index and donor processor (for flow -> donor = structure) ---*/
    /*--- Then we will pass on to the structural side the index (fea point) to the appropriate processor ---*/
    nBuffer_DonorIndices = 2 * MaxLocalVertexFlow;
    nBuffer_SetIndex = MaxLocalVertexStruct;
    
    /*--- Send and Recv buffers ---*/
    
    /*--- Buffers to send and receive the structural coordinates ---*/
    su2double *Buffer_Send_FlowTraction = new su2double[nBuffer_FlowTraction];
    su2double *Buffer_Recv_FlowTraction = NULL;
    
    /*--- Buffers to send and receive the donor index and processor ---*/
    long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
    long *Buffer_Recv_DonorIndices = NULL;
    
    /*--- Buffers to send and receive the new fluid coordinates ---*/
    su2double *Buffer_Send_StructTraction = NULL;
    su2double *Buffer_Recv_StructTraction = new su2double[nBuffer_StructTraction];
    
    /*--- Buffers to send and receive the fluid index ---*/
    long *Buffer_Send_SetIndex = NULL;
    long *Buffer_Recv_SetIndex = new long[nBuffer_SetIndex];
    
    /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/
    
    if (rank == MASTER_NODE) {
      Buffer_Recv_FlowTraction  = new su2double[size*nBuffer_FlowTraction];
      Buffer_Recv_DonorIndices = new long[size*nBuffer_DonorIndices];
      Buffer_Send_StructTraction = new su2double[size*nBuffer_StructTraction];
      Buffer_Send_SetIndex     = new long[size*nBuffer_SetIndex];
    }
    
    /*--- On the fluid side ---*/
    
    /*--- If this processor owns the marker we are looping at on the structural side ---*/
    
    /*--- First we initialize all of the indices and processors to -1 ---*/
    /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
      Buffer_Send_DonorIndices[iVertex] = -1;
    
    if (Marker_Flow >= 0) {
      
      /*--- We have identified the local index of the FEA marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < nLocalVertexFlow; iVertex++) {
        
        Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();
        
        Point_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorPoint();
        
        Processor_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorProcessor();
        
        // Get the normal at the vertex: this normal goes inside the fluid domain.
        Normal_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNormal();
        
        // Retrieve the values of pressure, viscosity and density
        if (incompressible) {
          
          Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
          Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
          
          if (viscous_flow) {
            
            Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
            Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
          }
        }
        else if (compressible) {
          
          Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
          Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
          
          if (viscous_flow) {
            
            Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
            Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
          }
        }
        
        // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
        for (iDim = 0; iDim < nDim; iDim++) {
          tn_f[iDim] = -(Pn-Pinf)*Normal_Flow[iDim];
        }
        
        // Calculate tn in the fluid nodes for the viscous term
        
        if ((incompressible || compressible) && viscous_flow) {
          
          // Divergence of the velocity
          div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
          if (incompressible) div_vel = 0.0;
          
          for (iDim = 0; iDim < nDim; iDim++) {
            
            for (jDim = 0 ; jDim < nDim; jDim++) {
              // Dirac delta
              Dij = 0.0; if (iDim == jDim) Dij = 1.0;
              
              // Viscous stress
              Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
              TWO3*Viscosity*div_vel*Dij;
              
              // Viscous component in the tn vector --> Units of force (non-dimensional).
              tn_f[iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
            }
          }
        }
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Buffer_Send_FlowTraction[iVertex*nDim+iDim] = tn_f[iDim]*factorForces*ModAmpl;
        }
        /*--- If this processor owns the node ---*/
        if (flow_geometry[MESH_0]->node[Point_Flow]->GetDomain()) {
          Buffer_Send_DonorIndices[2*iVertex]     = Point_Struct;
          Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Struct;
        }
        else {
          /*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
          Buffer_Send_DonorIndices[2*iVertex]     = -1;
          Buffer_Send_DonorIndices[2*iVertex + 1] = -1;
        }
        
      }
    }
    
    /*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
    SU2_MPI::Gather(Buffer_Send_FlowTraction, nBuffer_FlowTraction, MPI_DOUBLE, Buffer_Recv_FlowTraction, nBuffer_FlowTraction, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
    //    if (rank == MASTER_NODE) {
    //      cout << endl << "-----------------------------------------------------------" << endl;
    //      cout << "For tag " << iMarkerFSI << ":" << endl;
    //      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    //        cout << "The processor " << iProcessor << " has " << Buffer_Recv_nVertexStruct[iProcessor] << " nodes on the structural side and ";
    //        cout << Buffer_Recv_nVertexFlow[iProcessor] << " nodes on the fluid side " << endl;
    //      }
    //      cout << "The max number of vertices is " << MaxLocalVertexStruct << " on the structural side and ";
    //      cout << MaxLocalVertexFlow << " on the fluid side." << endl;
    //
    //      cout << "---------------- Check received buffers ---------------------" << endl;
    //      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    //        long initialIndex, initialIndex2;
    //        initialIndex = iProcessor*nBuffer_FlowTraction;
    //        initialIndex2 = iProcessor*nBuffer_DonorIndices;
    //        for (long iCheck = 0; iCheck < Buffer_Recv_nVertexStruct[iProcessor]; iCheck++) {
    //          cout << "From processor " << iProcessor << " we get coordinates (";
    //            for (iDim = 0; iDim < nDim; iDim++)
    //              cout << Buffer_Recv_FlowTraction[initialIndex+iCheck*nDim+iDim] << ",";
    //          cout << "), the donor index for the flow " << Buffer_Recv_DonorIndices[initialIndex2+iCheck*2] ;
    //          cout << " and the donor processor " << Buffer_Recv_DonorIndices[initialIndex2+iCheck*2+1] << endl;
    //
    //        }
    //      }
    //
    //    }
    
    /*--- Counter to determine where in the array we have to set the information ---*/
    long *Counter_Processor_Struct = NULL;
    long iProcessor_Flow = 0, iIndex_Flow = 0;
    long iProcessor_Struct = 0, iPoint_Struct = 0, iIndex_Struct = 0;
    long Point_Struct_Send, Processor_Struct_Send;
    
    /*--- Now we pack the information to send it over to the different processors ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- We set the counter to 0 ---*/
      Counter_Processor_Struct = new long[nProcessor];
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        Counter_Processor_Struct[iProcessor] = 0;
      }
      
      /*--- First we initialize the index vector to -1 ---*/
      /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
      for (iVertex = 0; iVertex < nProcessor*nBuffer_SetIndex; iVertex++)
        Buffer_Send_SetIndex[iVertex] = -2;
      
      /*--- As of now we do the loop over the flow points ---*/
      /*--- The number of points for flow and structure does not necessarily have to match ---*/
      /*--- In fact, it's possible that a processor asks for nStruct nodes and there are only ---*/
      /*--- nFlow < nStruct available; this is due to halo nodes ---*/
      
      /*--- For every processor from which we have received information ---*/
      /*--- (This is, for every processor on the structural side) ---*/
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        
        /*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
        iProcessor_Flow = iProcessor*nBuffer_FlowTraction;
        /*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
        iIndex_Flow = iProcessor*nBuffer_DonorIndices;
        
        /*--- For every vertex in the information retreived from iProcessor ---*/
        for (iVertex = 0; iVertex < Buffer_Recv_nVertexFlow[iProcessor]; iVertex++) {
          
          /*--- The processor and index for the flow are: ---*/
          Processor_Struct_Send = Buffer_Recv_DonorIndices[iIndex_Flow+iVertex*2+1];
          Point_Struct_Send     = Buffer_Recv_DonorIndices[iIndex_Flow+iVertex*2];
          
          /*--- Load the buffer at the appropriate position ---*/
          /*--- This is determined on the fluid side by:
           *--- Processor_Flow*nBuffer_StructTraction -> Initial position of the processor array (fluid side)
           *--- +
           *--- Counter_Processor_Struct*nDim -> Initial position of the nDim array for the particular point on the fluid side
           *--- +
           *--- iDim -> Position within the nDim array that corresponds to a point
           *---
           *--- While on the structural side is:
           *--- iProcessor*nBuffer_FlowTraction -> Initial position on the processor array (structural side)
           *--- +
           *--- iVertex*nDim -> Initial position of the nDim array for the particular point on the structural side
           */
          
          /*--- We check that we are not setting the value for a halo node ---*/
          if (Point_Struct_Send != -1) {
            iProcessor_Struct = Processor_Struct_Send*nBuffer_StructTraction;
            iIndex_Struct = Processor_Struct_Send*nBuffer_SetIndex;
            iPoint_Struct = Counter_Processor_Struct[Processor_Struct_Send]*nDim;
            
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_StructTraction[iProcessor_Struct + iPoint_Struct + iDim] = Buffer_Recv_FlowTraction[iProcessor_Flow + iVertex*nDim + iDim];
            
            /*--- We set the fluid index at an appropriate position matching the coordinates ---*/
            Buffer_Send_SetIndex[iIndex_Struct + Counter_Processor_Struct[Processor_Struct_Send]] = Point_Struct_Send;
            
            Counter_Processor_Struct[Processor_Struct_Send]++;
          }
          
        }
        
      }
      
      //      cout << "---------------- Check send buffers ---------------------" << endl;
      //
      //      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      //        long initialIndex, initialIndex2;
      //        initialIndex = iProcessor*nBuffer_StructTraction;
      //        initialIndex2 = iProcessor*nBuffer_SetIndex;
      //        for (long iCheck = 0; iCheck < Buffer_Recv_nVertexFlow[iProcessor]; iCheck++) {
      //          cout << "Processor " << iProcessor << " will receive the node " ;
      //          cout << Buffer_Send_SetIndex[initialIndex2+iCheck] << " which corresponds to the coordinates ";
      //          for (iDim = 0; iDim < nDim; iDim++)
      //            cout << "x" << iDim << "=" << Buffer_Send_StructTraction[initialIndex + iCheck*nDim + iDim] << ", ";
      //          cout << endl;
      //        }
      //
      //      }
      
    }
    
    /*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
    SU2_MPI::Scatter(Buffer_Send_StructTraction, nBuffer_StructTraction, MPI_DOUBLE, Buffer_Recv_StructTraction, nBuffer_StructTraction, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Scatter(Buffer_Send_SetIndex, nBuffer_SetIndex, MPI_LONG, Buffer_Recv_SetIndex, nBuffer_SetIndex, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
    long indexPoint_iVertex, Point_Struct_Check;
    long Point_Struct_Recv;
    
    /*--- For the flow marker we are studying ---*/
    if (Marker_Struct >= 0) {
      
      /*--- We have identified the local index of the Structural marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < nLocalVertexStruct; iVertex++) {
        
        Point_Struct_Recv = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetNode();
        
        if (fea_geometry[MESH_0]->node[Point_Struct_Recv]->GetDomain()) {
          /*--- Find the index of the point Point_Struct in the buffer Buffer_Recv_SetIndex ---*/
          indexPoint_iVertex = std::distance(Buffer_Recv_SetIndex, std::find(Buffer_Recv_SetIndex, Buffer_Recv_SetIndex + MaxLocalVertexStruct, Point_Struct_Recv));
          
          Point_Struct_Check = Buffer_Recv_SetIndex[indexPoint_iVertex];
          
          if (Point_Struct_Check < 0) {
            cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
            exit(EXIT_FAILURE);
          }
          
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[iDim] = Buffer_Recv_StructTraction[indexPoint_iVertex*nDim+iDim];
          
          /*--- Add to the Flow traction ---*/
          node[Point_Struct_Recv]->Add_FlowTraction(Residual);
          
        }
        
      }
      
    }
    
    delete [] Buffer_Send_FlowTraction;
    delete [] Buffer_Send_DonorIndices;
    delete [] Buffer_Recv_StructTraction;
    delete [] Buffer_Recv_SetIndex;
    
    if (rank == MASTER_NODE) {
      delete [] Buffer_Recv_nVertexStruct;
      delete [] Buffer_Recv_nVertexFlow;
      delete [] Buffer_Recv_FlowTraction;
      delete [] Buffer_Recv_DonorIndices;
      delete [] Buffer_Send_StructTraction;
      delete [] Buffer_Send_SetIndex;
      delete [] Counter_Processor_Struct;
    }
    
  }
  
#endif
  
  delete[] tn_f;
  
  
}

void CFEM_ElasticitySolver::SetFEA_Load_Int(CSolver ***flow_solution, CGeometry **fea_geometry,
                                            CGeometry **flow_geometry, CConfig *fea_config,
                                            CConfig *flow_config, CNumerics *fea_numerics) { }

void CFEM_ElasticitySolver::PredictStruct_Displacement(CGeometry **fea_geometry,
                                                       CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned short predOrder = fea_config->GetPredictorOrder();
  su2double Delta_t = fea_config->GetDelta_DynTime();
  unsigned long iPoint, iDim;
  su2double *solDisp, *solVel, *solVel_tn, *valPred;
  
  //To nPointDomain: we need to communicate the predicted solution after setting it
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    if (predOrder==0) fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    else if (predOrder==1) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim < nDim; iDim++) {
        valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
      }
      
    }
    else if (predOrder==2) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      solVel_tn = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim < nDim; iDim++) {
        valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
      }
      
    }
    else {
      cout<< "Higher order predictor not implemented. Solving with order 0." << endl;
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    }
  }
  
}

void CFEM_ElasticitySolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config,
                                                      CSolver ***fea_solution, unsigned long iFSIIter) {
  
  unsigned long iPoint, iDim;
  su2double rbuf_numAitk = 0, sbuf_numAitk = 0;
  su2double rbuf_denAitk = 0, sbuf_denAitk = 0;
  
  su2double *dispPred, *dispCalc, *dispPred_Old, *dispCalc_Old;
  su2double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
  su2double delta_deltaU[3] = {0.0, 0.0, 0.0};
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double Static_Time=fea_config->GetStatic_Time();
  su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn_Min, WAitkDyn;
  
  unsigned short RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  ofstream historyFile_FSI;
  bool writeHistFSI = fea_config->GetWrite_Conv_FSI();
  if (writeHistFSI && (rank == MASTER_NODE)) {
    char cstrFSI[200];
    string filenameHistFSI = fea_config->GetConv_FileName_FSI();
    strcpy (cstrFSI, filenameHistFSI.data());
    historyFile_FSI.open (cstrFSI, std::ios_base::app);
  }
  
  
  /*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/
  
  if (CurrentTime > Static_Time) {
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      
      if (writeHistFSI && (rank == MASTER_NODE)) {
        
        SetWAitken_Dyn(1.0);
        
        if (iFSIIter == 0) historyFile_FSI << " " << endl ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
        if (iFSIIter == 0) historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << 1.0 ;
        else historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << 1.0 << "," ;
      }
      
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      
      if (writeHistFSI && (rank == MASTER_NODE)) {
        
        SetWAitken_Dyn(fea_config->GetAitkenStatRelax());
        
        if (iFSIIter == 0) historyFile_FSI << " " << endl ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
        if (iFSIIter == 0) historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << fea_config->GetAitkenStatRelax() ;
        else historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << fea_config->GetAitkenStatRelax() << "," ;
      }
      
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      
      if (iFSIIter == 0) {
        
        WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
        WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();
        WAitkDyn_Min = fea_config->GetAitkenDynMinInit();
        
        WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);
        WAitkDyn = max(WAitkDyn, WAitkDyn_Min);
        
        SetWAitken_Dyn(WAitkDyn);
        if (writeHistFSI && (rank == MASTER_NODE)) {
          if (iFSIIter == 0) historyFile_FSI << " " << endl ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
          historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn ;
        }
        
      }
      else {
        // To nPointDomain; we need to communicate the values
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          
          dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
          dispPred_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred_Old();
          dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
          dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Old();
          
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*--- Compute the deltaU and deltaU_n+1 ---*/
            deltaU[iDim] = dispCalc_Old[iDim] - dispPred_Old[iDim];
            deltaU_p1[iDim] = dispCalc[iDim] - dispPred[iDim];
            
            /*--- Compute the difference ---*/
            delta_deltaU[iDim] = deltaU_p1[iDim] - deltaU[iDim];
            
            /*--- Add numerator and denominator ---*/
            sbuf_numAitk += deltaU[iDim] * delta_deltaU[iDim];
            sbuf_denAitk += delta_deltaU[iDim] * delta_deltaU[iDim];
            
          }
          
        }
        
#ifdef HAVE_MPI
        SU2_MPI::Allreduce(&sbuf_numAitk, &rbuf_numAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&sbuf_denAitk, &rbuf_denAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        rbuf_numAitk = sbuf_numAitk;
        rbuf_denAitk = sbuf_denAitk;
#endif
        
        WAitkDyn = GetWAitken_Dyn();
        
        if (rbuf_denAitk > 1E-15) {
          WAitkDyn = - 1.0 * WAitkDyn * rbuf_numAitk / rbuf_denAitk ;
        }
        
        WAitkDyn = max(WAitkDyn, 0.1);
        WAitkDyn = min(WAitkDyn, 1.0);
        
        SetWAitken_Dyn(WAitkDyn);
        
        if (writeHistFSI && (rank == MASTER_NODE)) {
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
          historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn << "," ;
        }
        
      }
      
    }
    else {
      if (rank == MASTER_NODE) cout << "No relaxation method used. " << endl;
    }
    
  }
  
  if (writeHistFSI && (rank == MASTER_NODE)) {historyFile_FSI.close();}
  
}

void CFEM_ElasticitySolver::SetAitken_Relaxation(CGeometry **fea_geometry,
                                                 CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint, iDim;
  unsigned short RelaxMethod_FSI;
  su2double *dispPred, *dispCalc;
  su2double WAitken;
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double Static_Time=fea_config->GetStatic_Time();
  
  RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  /*--- Only when there is movement it makes sense to update the solutions... ---*/
  
  if (CurrentTime > Static_Time) {
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      WAitken = 1.0;
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      WAitken = fea_config->GetAitkenStatRelax();
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      WAitken = GetWAitken_Dyn();
    }
    else {
      WAitken = 1.0;
    }
    
    // To nPointDomain; we need to communicate the solutions (predicted, old and old predicted) after this routine
    for (iPoint=0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve pointers to the predicted and calculated solutions ---*/
      dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      
      /*--- Set predicted solution as the old predicted solution ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred_Old();
      
      /*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Old(dispCalc);
      
      /*--- Apply the Aitken relaxation ---*/
      for (iDim=0; iDim < nDim; iDim++) {
        dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Update_StructSolution(CGeometry **fea_geometry,
                                                  CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint;
  su2double *valSolutionPred;
  
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    
    valSolutionPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
    
    fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution(valSolutionPred);
    
  }
  
  /*--- Perform the MPI communication of the solution, displacements only ---*/
  
  Set_MPI_Solution_DispOnly(fea_geometry[MESH_0], fea_config);
  
}


void CFEM_ElasticitySolver::Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  unsigned long nTotalPoint = 1;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

  unsigned long ExtIter = config->GetExtIter();

  su2double reference_geometry = 0.0, current_solution = 0.0;

  bool fsi = config->GetFSI_Simulation();

  su2double objective_function = 0.0, objective_function_reduce = 0.0;
  su2double weight_OF = 1.0;

  su2double objective_function_averaged = 0.0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nPointDomain,  &nTotalPoint,  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nTotalPoint        = nPointDomain;
#endif

  weight_OF = config->GetRefGeom_Penalty() / nTotalPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    for (iVar = 0; iVar < nVar; iVar++){

      /*--- Retrieve the value of the reference geometry ---*/
      reference_geometry = node[iPoint]->GetReference_Geometry(iVar);

      /*--- Retrieve the value of the current solution ---*/
      current_solution = node[iPoint]->GetSolution(iVar);

      /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
      objective_function += weight_OF * (current_solution - reference_geometry)*(current_solution - reference_geometry);
    }

  }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&objective_function,  &objective_function_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    objective_function_reduce        = objective_function;
#endif

  Total_OFRefGeom = objective_function_reduce + PenaltyValue;

  Global_OFRefGeom += Total_OFRefGeom;
  objective_function_averaged = Global_OFRefGeom / (ExtIter + 1.0 + EPS);

  bool direct_diff = ((config->GetDirectDiff() == D_YOUNG) ||
                      (config->GetDirectDiff() == D_POISSON) ||
                      (config->GetDirectDiff() == D_RHO) ||
                      (config->GetDirectDiff() == D_RHO_DL) ||
                      (config->GetDirectDiff() == D_EFIELD) ||
                      (config->GetDirectDiff() == D_MACH) ||
                      (config->GetDirectDiff() == D_PRESSURE));

  if ((direct_diff) && (rank == MASTER_NODE)){

    ofstream myfile_res;

    if (config->GetDirectDiff() == D_YOUNG) myfile_res.open ("Output_Direct_Diff_E.txt", ios::app);
    if (config->GetDirectDiff() == D_POISSON) myfile_res.open ("Output_Direct_Diff_Nu.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO) myfile_res.open ("Output_Direct_Diff_Rho.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO_DL) myfile_res.open ("Output_Direct_Diff_Rho_DL.txt", ios::app);
    if (config->GetDirectDiff() == D_EFIELD) myfile_res.open ("Output_Direct_Diff_EField.txt", ios::app);

    if (config->GetDirectDiff() == D_MACH) myfile_res.open ("Output_Direct_Diff_Mach.txt");
    if (config->GetDirectDiff() == D_PRESSURE) myfile_res.open ("Output_Direct_Diff_Pressure.txt");

    myfile_res.precision(15);

    myfile_res << scientific << Total_OFRefGeom << "\t";

    if (dynamic) myfile_res << scientific << objective_function_averaged << "\t";

    su2double local_forward_gradient = 0.0;
    su2double averaged_gradient = 0.0;

    local_forward_gradient = SU2_TYPE::GetDerivative(Total_OFRefGeom);

    if (fsi) {
      Total_ForwardGradient = local_forward_gradient;
      averaged_gradient     = Total_ForwardGradient / (ExtIter + 1.0);
    }
    else {
      Total_ForwardGradient += local_forward_gradient;
      averaged_gradient      = Total_ForwardGradient / (ExtIter + 1.0);
    }

    myfile_res << scientific << local_forward_gradient << "\t";

    myfile_res << scientific << averaged_gradient << "\t";

    myfile_res << endl;

    myfile_res.close();

    if (config->GetDirectDiff() == D_YOUNG)   cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the Young Modulus: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_POISSON) cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the Poisson's ratio: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO)     cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the structural density: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO_DL)  cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the dead weight: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_EFIELD)  cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the electric field: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_MACH)    cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the Mach number: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_PRESSURE) cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the freestream Pressure: " << Total_ForwardGradient << "." << endl;
  }
  else
  {

      // TODO: Temporary output file for the objective function. Will be integrated in the output once is refurbished.
   if (rank == MASTER_NODE){
      cout << "Objective function: " << Total_OFRefGeom << "." << endl;
      ofstream myfile_res;
      myfile_res.open ("of_refgeom.dat");
      myfile_res.precision(15);
      if (dynamic) myfile_res << scientific << objective_function_averaged << endl;
      else myfile_res << scientific << Total_OFRefGeom << endl;
      myfile_res.close();
      if (fsi){
          ofstream myfile_his;
          myfile_his.open ("history_refgeom.dat",ios::app);
          myfile_his.precision(15);
          myfile_his << scientific << Total_OFRefGeom << endl;
          myfile_his.close();
      }
    }

  }

}

void CFEM_ElasticitySolver::Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

  unsigned long ExtIter = config->GetExtIter();

  su2double reference_geometry = 0.0, current_solution = 0.0;

  bool fsi = config->GetFSI_Simulation();

  su2double objective_function = 0.0, objective_function_reduce = 0.0;
  su2double distance_sq = 0.0 ;
  su2double weight_OF = 1.0;

  su2double objective_function_averaged = 0.0;

  /*--- TEMPORARY, for application in dynamic TestCase ---*/
  su2double difX = 0.0, difX_reduce = 0.0;
  su2double difY = 0.0, difY_reduce = 0.0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  weight_OF = config->GetRefNode_Penalty();

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    if (geometry->node[iPoint]->GetGlobalIndex() == config->GetRefNode_ID() ){

      for (iVar = 0; iVar < nVar; iVar++){

        /*--- Retrieve the value of the reference geometry ---*/
        reference_geometry = config->GetRefNode_Displacement(iVar);

        /*--- Retrieve the value of the current solution ---*/
        current_solution = node[iPoint]->GetSolution(iVar);

        /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
        distance_sq +=  (current_solution - reference_geometry)*(current_solution - reference_geometry);
      }

      objective_function = weight_OF * sqrt(distance_sq);

      difX = node[iPoint]->GetSolution(0) - config->GetRefNode_Displacement(0);
      difY = node[iPoint]->GetSolution(1) - config->GetRefNode_Displacement(1);

    }

  }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&objective_function,  &objective_function_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&difX,  &difX_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&difY,  &difY_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    objective_function_reduce        = objective_function;
    difX_reduce                      = difX;
    difY_reduce                      = difY;
#endif

  Total_OFRefNode = objective_function_reduce + PenaltyValue;

  Global_OFRefNode += Total_OFRefNode;
  objective_function_averaged = Global_OFRefNode / (ExtIter + 1.0 + EPS);

  bool direct_diff = ((config->GetDirectDiff() == D_YOUNG) ||
                      (config->GetDirectDiff() == D_POISSON) ||
                      (config->GetDirectDiff() == D_RHO) ||
                      (config->GetDirectDiff() == D_RHO_DL) ||
                      (config->GetDirectDiff() == D_EFIELD) ||
                      (config->GetDirectDiff() == D_MACH) ||
                      (config->GetDirectDiff() == D_PRESSURE));

  if ((direct_diff) && (rank == MASTER_NODE)){

    ofstream myfile_res;

    if (config->GetDirectDiff() == D_YOUNG) myfile_res.open ("Output_Direct_Diff_E.txt", ios::app);
    if (config->GetDirectDiff() == D_POISSON) myfile_res.open ("Output_Direct_Diff_Nu.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO) myfile_res.open ("Output_Direct_Diff_Rho.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO_DL) myfile_res.open ("Output_Direct_Diff_Rho_DL.txt", ios::app);
    if (config->GetDirectDiff() == D_EFIELD) myfile_res.open ("Output_Direct_Diff_EField.txt", ios::app);

    if (config->GetDirectDiff() == D_MACH) myfile_res.open ("Output_Direct_Diff_Mach.txt");
    if (config->GetDirectDiff() == D_PRESSURE) myfile_res.open ("Output_Direct_Diff_Pressure.txt");

    myfile_res.precision(15);

    myfile_res << scientific << Total_OFRefNode << "\t";

    if (dynamic) myfile_res << scientific << objective_function_averaged << "\t";

    su2double local_forward_gradient = 0.0;
    su2double averaged_gradient = 0.0;

    local_forward_gradient = SU2_TYPE::GetDerivative(Total_OFRefNode);

    if (fsi) {
      Total_ForwardGradient = local_forward_gradient;
      averaged_gradient     = Total_ForwardGradient / (ExtIter + 1.0);
    }
    else {
      Total_ForwardGradient += local_forward_gradient;
      averaged_gradient      = Total_ForwardGradient / (ExtIter + 1.0);
    }

    myfile_res << scientific << local_forward_gradient << "\t";

    myfile_res << scientific << averaged_gradient << "\t";

    myfile_res << scientific << difX_reduce << "\t";
    myfile_res << scientific << difY_reduce << "\t";

    myfile_res << endl;

    myfile_res.close();

    if (config->GetDirectDiff() == D_YOUNG)   cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the Young Modulus: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_POISSON) cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the Poisson's ratio: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO)     cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the structural density: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO_DL)  cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the dead weight: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_EFIELD)  cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the electric field: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_MACH)    cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the Mach number: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_PRESSURE) cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the freestream Pressure: " << Total_ForwardGradient << "." << endl;
  }
  else
  {

    // TODO: Temporary output file for the objective function. Will be integrated in the output once is refurbished.
    if (rank == MASTER_NODE){
      cout << "Objective function: " << Total_OFRefNode << "." << endl;
      ofstream myfile_res;
      myfile_res.open ("of_refnode.dat");
      myfile_res.precision(15);
      if (dynamic) myfile_res << scientific << objective_function_averaged << endl;
      else myfile_res << scientific << Total_OFRefNode << endl;
      myfile_res.close();

      ofstream myfile_his;
      myfile_his.open ("history_refnode.dat",ios::app);
      myfile_his.precision(15);
      myfile_his << ExtIter << "\t";
      myfile_his << scientific << Total_OFRefNode << "\t";
      myfile_his << scientific << objective_function_averaged << "\t";
      myfile_his << scientific << difX_reduce << "\t";
      myfile_his << scientific << difY_reduce << endl;
      myfile_his.close();

    }

  }


}

void CFEM_ElasticitySolver::Stiffness_Penalty(CGeometry *geometry, CSolver **solver, CNumerics **numerics, CConfig *config){

  unsigned long iElem;
  unsigned short iNode, iDim, nNodes = 0;

  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol;
  int EL_KIND = 0;

  su2double elementVolume, dvValue, ratio;
  su2double weightedValue = 0.0;
  su2double weightedValue_reduce = 0.0;
  su2double totalVolume = 0.0;
  su2double totalVolume_reduce = 0.0;

  /*--- Loops over the elements in the domain ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}

    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; iDim++) {
            val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
            val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
            element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
            element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        }
    }

    // Avoid double-counting elements:
    // Only add the value if the first node is in the domain
    if (geometry->node[indexNode[0]]->GetDomain()){

        // Compute the area/volume of the element
        if (nDim == 2)
        	elementVolume = element_container[FEA_TERM][EL_KIND]->ComputeArea();
        else
        	elementVolume = element_container[FEA_TERM][EL_KIND]->ComputeVolume();

        // Compute the total volume
        totalVolume += elementVolume;

        // Retrieve the value of the design variable
        dvValue = numerics[FEA_TERM]->Get_DV_Val(element_properties[iElem]->GetDV());

        // Add the weighted sum of the value of the design variable
        weightedValue += dvValue * elementVolume;

    }

  }

// Reduce value across processors for parallelization

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&weightedValue,  &weightedValue_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&totalVolume,  &totalVolume_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    weightedValue_reduce        = weightedValue;
    totalVolume_reduce          = totalVolume;
#endif

    ratio = 1.0 - weightedValue_reduce/totalVolume_reduce;

    PenaltyValue = config->GetTotalDV_Penalty() * ratio * ratio;



}

void CFEM_ElasticitySolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  unsigned short iVar, nSolVar;
  unsigned long index;

  ifstream restart_file;
  string restart_filename, filename, text_line;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry[MESH_0]->GetnZone();

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fluid_structure = config->GetFSI_Simulation();
  bool discrete_adjoint = config->GetDiscrete_Adjoint();

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

  su2double *Sol = new su2double[nSolVar];

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Restart the solution from file information ---*/

  filename = config->GetSolution_FEMFileName();

  /*--- If multizone, append zone name ---*/

  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  if (dynamic) {
    filename = config->GetUnsteady_FileName(filename, val_iter);
  }

  /*--- Read all lines in the restart file ---*/

  int counter = 0;
  long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, filename);
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
      for (iVar = 0; iVar < nSolVar; iVar++) Sol[iVar] = Restart_Data[index+iVar];

      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint_Local]->SetSolution(iVar, Sol[iVar]);
        if (dynamic) {
          node[iPoint_Local]->SetSolution_time_n(iVar, Sol[iVar]);
          node[iPoint_Local]->SetSolution_Vel(iVar, Sol[iVar+nVar]);
          node[iPoint_Local]->SetSolution_Vel_time_n(iVar, Sol[iVar+nVar]);
          node[iPoint_Local]->SetSolution_Accel(iVar, Sol[iVar+2*nVar]);
          node[iPoint_Local]->SetSolution_Accel_time_n(iVar, Sol[iVar+2*nVar]);
        }
        if (fluid_structure && !dynamic) {
          node[iPoint_Local]->SetSolution_Pred(iVar, Sol[iVar]);
          node[iPoint_Local]->SetSolution_Pred_Old(iVar, Sol[iVar]);
        }
        if (fluid_structure && discrete_adjoint){
          node[iPoint_Local]->SetSolution_Old(iVar, Sol[iVar]);
        }
      }
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
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

  /*--- MPI. If dynamic, we also need to communicate the old solution ---*/

  solver[MESH_0][FEA_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  if (dynamic) solver[MESH_0][FEA_SOL]->Set_MPI_Solution_Old(geometry[MESH_0], config);
  if (fluid_structure && !dynamic){
      solver[MESH_0][FEA_SOL]->Set_MPI_Solution_Pred(geometry[MESH_0], config);
      solver[MESH_0][FEA_SOL]->Set_MPI_Solution_Pred_Old(geometry[MESH_0], config);
  }

  delete [] Sol;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;
  
}
