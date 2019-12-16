/*!
 * \file solver_direct_elasticity.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
 * \author R. Sanchez
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
#include "../include/variables/CFEABoundVariable.hpp"
#include "../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../Common/include/omp_structure.hpp"
#include <algorithm>

CFEASolver::CFEASolver(bool mesh_deform_mode) : CSolver(mesh_deform_mode) {

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
  unsigned short iTerm;
  element_container = new CElement** [MAX_TERMS]();
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS*omp_get_max_threads()]();

  nodes = nullptr;

  element_properties = NULL;

  SolRest = NULL;

  iElem_iDe = NULL;

  topol_filter_applied = false;

  element_based = false;
}

CFEASolver::CFEASolver(CGeometry *geometry, CConfig *config) : CSolver() {

  unsigned long iPoint;
  unsigned short iVar;
  unsigned short iTerm;

  bool dynamic = (config->GetTime_Domain());

  /*--- Test whether we consider dielectric elastomers ---*/
  bool de_effects = config->GetDE_Effects();

  /*--- A priori we don't have an element-based input file (most of the applications will be like this) ---*/
  element_based = false;

  topol_filter_applied = false;

  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();

  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();

  /*--- Number of different terms for FEA ---*/
  nFEA_Terms = 1;
  if (de_effects) nFEA_Terms++; // The DE term is DE_TERM = 1

  /*--- Here is where we assign the kind of each element ---*/

  /*--- First level: different possible terms of the equations ---*/
  element_container = new CElement** [MAX_TERMS]();
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS*omp_get_max_threads()]();

  if (nDim == 2) {
    for(int thread = 0; thread < omp_get_max_threads(); ++thread) {

      const int offset = thread*MAX_FE_KINDS;

      /*--- Basic terms ---*/
      element_container[FEA_TERM][EL_TRIA+offset] = new CTRIA1();
      element_container[FEA_TERM][EL_QUAD+offset] = new CQUAD4();

      if (de_effects){
        element_container[DE_TERM][EL_TRIA+offset] = new CTRIA1();
        element_container[DE_TERM][EL_QUAD+offset] = new CQUAD4();
      }
    }
  }
  else {
    for(int thread = 0; thread < omp_get_max_threads(); ++thread) {

      const int offset = thread*MAX_FE_KINDS;

      element_container[FEA_TERM][EL_TETRA+offset] = new CTETRA1();
      element_container[FEA_TERM][EL_HEXA +offset] = new CHEXA8 ();
      element_container[FEA_TERM][EL_PYRAM+offset] = new CPYRAM5();
      element_container[FEA_TERM][EL_PRISM+offset] = new CPRISM6();

      if (de_effects){
        element_container[DE_TERM][EL_TETRA+offset] = new CTETRA1();
        element_container[DE_TERM][EL_HEXA +offset] = new CHEXA8 ();
        element_container[DE_TERM][EL_PYRAM+offset] = new CPYRAM5();
        element_container[DE_TERM][EL_PRISM+offset] = new CPRISM6();
      }
    }
  }

  /*--- Set element properties ---*/
  Set_ElementProperties(geometry, config);

  Total_CFEA        = 0.0;
  WAitken_Dyn       = 0.0;
  WAitken_Dyn_tn1   = 0.0;
  loadIncrement     = 0.0;

  SetFSI_ConvValue(0,0.0);
  SetFSI_ConvValue(1,0.0);

  nVar = nDim;

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual = new su2double[nVar]();
  Residual_RMS = new su2double[nVar]();
  Residual_Max = new su2double[nVar]();
  Point_Max = new unsigned long[nVar]();
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim]();
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution = new su2double[nVar]();

  /*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/

  unsigned short nSolVar;
  string text_line, filename;
  ifstream restart_file;

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

  SolRest = new su2double[nSolVar]();

  /*--- Initialize from zero everywhere ---*/

  nodes = new CFEABoundVariable(SolRest, nPoint, nDim, nVar, config);
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
  static_cast<CFEABoundVariable*>(nodes)->AllocateBoundaryVariables(config);


  bool reference_geometry = config->GetRefGeom();
  if (reference_geometry) Set_ReferenceGeometry(geometry, config);

  bool prestretch_fem = config->GetPrestretch();
  if (prestretch_fem) Set_Prestretch(geometry, config);

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

#ifdef HAVE_OMP
  /*--- Get the element coloring. ---*/

  const auto& coloring = geometry->GetElementColoring();

  auto nColor = coloring.getOuterSize();
  ElemColoring.resize(nColor);

  for(auto iColor = 0ul; iColor < nColor; ++iColor) {
    ElemColoring[iColor].size = coloring.getNumNonZeros(iColor);
    ElemColoring[iColor].indices = coloring.innerIdx(iColor);
  }

  ColorGroupSize = 1; /// TODO: This needs to come from geometry or config

  omp_chunk_size = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_SIZE);
#endif

  iElem_iDe = NULL;

  /*--- Initialize the value of the total objective function ---*/
  Total_OFRefGeom = 0.0;
  Total_OFRefNode = 0.0;
  Total_OFVolFrac = 0.0;

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

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (config->GetMultizone_Residual()){

    FSI_Residual      = 0.0;
    RelaxCoeff        = 1.0;
    ForceCoeff        = 1.0;

    Residual_BGS      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 1.0;
    Residual_Max_BGS  = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar]();
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim]();
    }
  }
  else{
    ForceCoeff = 1.0;
  }

  /*--- Penalty value - to maintain constant the stiffness in optimization problems - TODO: this has to be improved ---*/
  PenaltyValue = 0.0;

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION_FEA);
  CompleteComms(geometry, config, SOLUTION_FEA);

  /*--- If dynamic, we also need to communicate the old solution ---*/

  if(dynamic) {
    InitiateComms(geometry, config, SOLUTION_FEA_OLD);
    CompleteComms(geometry, config, SOLUTION_FEA_OLD);
  }

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "FEA";
}

CFEASolver::~CFEASolver(void) {

  unsigned short iVar, jVar;
  unsigned long iElem;

  if (element_container != NULL) {
    for (iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (jVar = 0; jVar < MAX_FE_KINDS*omp_get_max_threads(); jVar++) {
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

  delete [] SolRest;

  if (iElem_iDe != NULL) delete [] iElem_iDe;

  if (nodes != nullptr) delete nodes;
}

void CFEASolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) {

  unsigned long iElem;
  unsigned long index;
  unsigned long elProperties[4];

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  bool topology_mode = config->GetTopology_Optimization();

  string filename;
  ifstream properties_file;

  element_properties = new CProperty*[nElement];

  /*--- Restart the solution from file information ---*/

  filename = config->GetFEA_FileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone, ".dat");

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {
    if (rank == MASTER_NODE){
      cout << "There is no element-based properties file." << endl;
      cout << "The structural domain has uniform properties." << endl;

      if (topology_mode)
        SU2_MPI::Error("Topology mode requires an element-based properties file.",CURRENT_FUNCTION);
    }

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

        /*--- For backwards compatibility we only read a fifth column in topology mode ---*/
        if (topology_mode) {
          su2double elDensity;
          point_line >> elDensity;
          element_properties[iElem_Local]->SetDesignDensity(elDensity);
        }

        iElem_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iElem_Global_Local < nElement) { sbuf_NotMatching = 1; }

    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

    if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The properties file ") + filename + string(" doesn't match with the mesh file!\n")  +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Close the restart file ---*/

    properties_file.close();

    /*--- Free memory needed for the transformation ---*/

    delete [] Global2Local;

  }

}


void CFEASolver::Set_Prestretch(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned long index;

  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  string filename;
  ifstream prestretch_file;


  /*--- Restart the solution from file information ---*/

  filename = config->GetPrestretch_FEMFileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone, ".dat");

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  prestretch_file.open(filename.data(), ios::in);

  /*--- In case there is no file ---*/

  if (prestretch_file.fail()) {
    SU2_MPI::Error(string("There is no FEM prestretch reference file ") + filename, CURRENT_FUNCTION);
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

      for (iVar = 0; iVar < nVar; iVar++) nodes->SetPrestretch(iPoint_Local,iVar, Solution[iVar]);

      iPoint_Global_Local++;
    }
    iPoint_Global++;
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Close the restart file ---*/

  prestretch_file.close();

#ifdef HAVE_MPI
  /*--- We need to communicate here the prestretched geometry for the halo nodes. ---*/
  /*--- We avoid creating a new function as this may be reformatted.              ---*/

  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

  int send_to, receive_from;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = nodes->GetPrestretch(iPoint,iVar);
      }

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, NULL);

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
          nodes->SetPrestretch(iPoint,iVar, SolRest[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }
#endif

}

void CFEASolver::Set_ReferenceGeometry(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;

  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short file_format = config->GetRefGeom_FileFormat();

  string filename;
  ifstream reference_file;


  /*--- Restart the solution from file information ---*/

  filename = config->GetRefGeom_FEMFileName();

  /*--- If multizone, append zone name ---*/

  filename = config->GetMultizone_FileName(filename, iZone, ".csv");

  reference_file.open(filename.data(), ios::in);

  /*--- In case there is no file ---*/

  if (reference_file.fail()) {
    SU2_MPI::Error( "There is no FEM reference geometry file!!", CURRENT_FUNCTION);
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

    vector<string> point_line = PrintingToolbox::split(text_line, ',');

    /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/

    iPoint_Local = Global2Local[iPoint_Global];

    if (iPoint_Local >= 0) {

      if (nDim == 2){
        Solution[0] = PrintingToolbox::stod(point_line[3]);
        Solution[1] = PrintingToolbox::stod(point_line[4]);
      } else {
        Solution[0] = PrintingToolbox::stod(point_line[4]);
        Solution[1] = PrintingToolbox::stod(point_line[5]);
        Solution[2] = PrintingToolbox::stod(point_line[6]);
      }

      for (iVar = 0; iVar < nVar; iVar++) nodes->SetReference_Geometry(iPoint_Local,iVar, Solution[iVar]);

      iPoint_Global_Local++;
    }
    iPoint_Global++;
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n")  +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- I don't think we need to communicate ---*/

  /*--- Close the restart file ---*/

  reference_file.close();

  /*--- Free memory needed for the transformation ---*/

  delete [] Global2Local;

}


void CFEASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics,
                               unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;
  bool initial_calc = (config->GetTimeIter() == 0) && (config->GetInnerIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetInnerIter() == 0);                          // Checks if it is the first iteration
  bool dynamic = (config->GetTime_Domain());              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool restart = config->GetRestart();                        // Restart analysis
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetTimeIter()) ==SU2_TYPE::Int(config->GetRestart_Iter())); // Initial calculation for restart

  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);     // Discrete adjoint FEM solver

  bool incremental_load = config->GetIncrementalLoad();               // If an incremental load is applied

  bool body_forces = config->GetDeadLoad();                     // Body forces (dead loads).

  bool fsi = config->GetFSI_Simulation();
  bool consistent_interpolation = (!config->GetConservativeInterpolation() ||
                                  (config->GetKindInterpolation() == WEIGHTED_AVERAGE));

  bool topology_mode = config->GetTopology_Optimization();  // Density-based topology optimization


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
   * For dynamic problems we also need to recompute the Jacobian as that is where the RHS is computed
   * as a residual (and we need that for AD).
   */
  if ((initial_calc && linear_analysis)||
      (restart && initial_calc_restart && linear_analysis) ||
      (dynamic && disc_adj_fem) ||
      (dynamic && linear_analysis)) {
    Jacobian.SetValZero();
  }

  /*
   * For topology optimization we apply a filter on the design density field to avoid
   * numerical issues (checkerboards), ensure mesh independence, and impose a length scale.
   * This has to be done before computing the mass matrix and the dead load terms.
   * This filter, and the volume fraction objective function, require the element volumes,
   * so we ask "geometry" to compute them.
   * This only needs to be done for the undeformed (initial) shape.
   */
  if (topology_mode && !topol_filter_applied) {
    geometry->SetElemVolume(config);
    FilterElementDensities(geometry,config);
    topol_filter_applied = true;
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
    Compute_MassMatrix(geometry, numerics, config);
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
      for (iPoint = 0; iPoint < nPoint; iPoint++) nodes->Clear_BodyForces_Res(iPoint);
    }
    // Compute the dead load term
    Compute_DeadLoad(geometry, numerics, config);
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
    for (iPoint = 0; iPoint < nPoint; iPoint++) nodes->Clear_SurfaceLoad_Res(iPoint);
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
            nodes->Clear_SurfaceLoad_Res(iPoint);
          }
        }
        break;
      case DAMPER_BOUNDARY:
        /*--- For all the vertices in the marker iMarker ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          /*--- Retrieve the point ID ---*/
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          /*--- Clear the residual of the node, to avoid adding on previous values ---*/
          nodes->Clear_SurfaceLoad_Res(iPoint);
        }
        break;
    }
  }

  /*
   * FSI loads (computed upstream) need to be integrated if a nonconservative interpolation scheme is in use
   */
  if (fsi && first_iter && consistent_interpolation) Integrate_FSI_Loads(geometry,config);

}

void CFEASolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEASolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/

  if (config->GetIncrementalLoad()) nodes->Set_OldSolution();

}

void CFEASolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/

  if (config->GetIncrementalLoad()) nodes->Set_Solution();

}

void CFEASolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, jNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {

          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

          for (iDim = 0; iDim < nDim; iDim++) {
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            su2double val_Sol = nodes->GetSolution(indexNode[iNode],iDim) + val_Coord;
            element->SetRef_Coord(iNode, iDim, val_Coord);
            element->SetCurr_Coord(iNode, iDim, val_Sol);
          }
        }

        /*--- In topology mode determine the penalty to apply to the stiffness. ---*/
        su2double simp_penalty = 1.0;
        if (topology_mode) {
          su2double density = element_properties[iElem]->GetPhysicalDensity();
          simp_penalty = simp_minstiff+(1.0-simp_minstiff)*pow(density,simp_exponent);
        }

        /*--- Set the properties of the element ---*/
        element->Set_ElProperties(element_properties[iElem]);

        /*--- Compute the components of the jacobian and the stress term, one numerics per thread. ---*/
        int NUM_TERM = thread*MAX_TERMS + element_based? element_properties[iElem]->GetMat_Mod() : FEA_TERM;

        numerics[NUM_TERM]->Compute_Tangent_Matrix(element, config);

        /*--- Update residual and stiffness matrix with contributions from the element. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {

          auto Ta = element->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes[indexNode[iNode]*nVar+iVar] -= simp_penalty*Ta[iVar];

          for (jNode = 0; jNode < nNodes; jNode++) {
            auto Kab = element->Get_Kab(iNode, jNode);
            Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], simp_penalty, Kab);
          }
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  const bool prestretch_fem = config->GetPrestretch();
  const bool de_effects = config->GetDE_Effects();

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, jNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* fea_elem = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];
        CElement* de_elem = element_container[DE_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {

          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

          for (iDim = 0; iDim < nDim; iDim++) {
            /*--- Compute current coordinate. ---*/
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            su2double val_Sol = nodes->GetSolution(indexNode[iNode],iDim) + val_Coord;

            /*--- If pre-stretched the reference coordinate is stored in the nodes. ---*/
            if (prestretch_fem)
              val_Coord = nodes->GetPrestretch(indexNode[iNode],iDim);

            /*--- Set coordinates. ---*/
            fea_elem->SetCurr_Coord(iNode, iDim, val_Sol);
            fea_elem->SetRef_Coord(iNode, iDim, val_Coord);

            if (de_effects) {
              de_elem->SetCurr_Coord(iNode, iDim, val_Sol);
              de_elem->SetRef_Coord(iNode, iDim, val_Coord);
            }
          }
        }

        /*--- In topology mode determine the penalty to apply to the stiffness. ---*/
        su2double simp_penalty = 1.0;
        if (topology_mode) {
          su2double density = element_properties[iElem]->GetPhysicalDensity();
          simp_penalty = simp_minstiff+(1.0-simp_minstiff)*pow(density,simp_exponent);
        }

        /*--- Set the properties of the element. ---*/
        fea_elem->Set_ElProperties(element_properties[iElem]);
        if (de_effects)
          de_elem->Set_ElProperties(element_properties[iElem]);

        /*--- Compute the components of the Jacobian and the stress term for the material. ---*/
        int NUM_TERM = thread*MAX_TERMS + element_based? element_properties[iElem]->GetMat_Mod() : FEA_TERM;

        numerics[NUM_TERM]->Compute_Tangent_Matrix(fea_elem, config);

        /*--- Compute the electric component of the Jacobian and the stress term. ---*/
        if (de_effects)
          numerics[DE_TERM + thread*MAX_TERMS]->Compute_Tangent_Matrix(de_elem, config);

        /*--- Update residual and stiffness matrix with contributions from the element. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {

          auto Ta = fea_elem->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes[indexNode[iNode]*nVar+iVar] -= simp_penalty*Ta[iVar];

          /*--- Retrieve the electric contribution to the residual. ---*/
          if (de_effects) {
            auto Ta_DE = de_elem->Get_Kt_a(iNode);
            for (iVar = 0; iVar < nVar; iVar++)
              LinSysRes[indexNode[iNode]*nVar+iVar] -= simp_penalty*Ta_DE[iVar];
          }

          for (jNode = 0; jNode < nNodes; jNode++) {

            /*--- Get a pointer to the matrix block to perform the update. ---*/
            auto Kij = Jacobian.GetBlock(indexNode[iNode], indexNode[jNode]);

            /*--- Retrieve the values of the FEA term. ---*/
            auto Kab = fea_elem->Get_Kab(iNode, jNode);
            su2double Ks_ab = fea_elem->Get_Ks_ab(iNode, jNode);

            /*--- Full block. ---*/
            for (iVar = 0; iVar < nVar*nVar; iVar++)
              Kij[iVar] += SU2_TYPE::GetValue(simp_penalty*Kab[iVar]);

            /*--- Only the block's diagonal. ---*/
            for (iVar = 0; iVar < nVar; iVar++)
              Kij[iVar*(nVar+1)] += SU2_TYPE::GetValue(simp_penalty*Ks_ab);

            /*--- Retrieve the electric contribution to the Jacobian ---*/
            if (de_effects) {
              //auto Kab_DE = de_elem->Get_Kab(iNode, jNode);
              su2double Ks_ab_DE = de_elem->Get_Ks_ab(iNode, jNode);

              for (iVar = 0; iVar < nVar; iVar++)
                Kij[iVar*(nVar+1)] += SU2_TYPE::GetValue(simp_penalty*Ks_ab_DE);
            }
          }
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_MassMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, jNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {
          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
          for (iDim = 0; iDim < nDim; iDim++) {
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            element->SetRef_Coord(iNode, iDim, val_Coord);
          }
        }

        /*--- In topology mode determine the penalty to apply to the mass,
         *    linear function of the physical density. ---*/
        su2double simp_penalty = 1.0;
        if (topology_mode) {
          simp_penalty = simp_minstiff+(1.0-simp_minstiff)*element_properties[iElem]->GetPhysicalDensity();
        }

        /*--- Set the properties of the element and compute its mass matrix. ---*/
        element->Set_ElProperties(element_properties[iElem]);

        numerics[FEA_TERM + thread*MAX_TERMS]->Compute_Mass_Matrix(element, config);

        /*--- Add contributions of this element to the mass matrix. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {

          for (jNode = 0; jNode < nNodes; jNode++) {

            auto Mij = MassMatrix.GetBlock(indexNode[iNode], indexNode[jNode]);
            su2double Mab = element->Get_Mab(iNode, jNode);

            for (iVar = 0; iVar < nVar; iVar++)
              Mij[iVar*(nVar+1)] += simp_penalty*Mab;
          }
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_MassRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    TimeRes.SetValZero();
    SU2_OMP_BARRIER

#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, jNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {
          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
          for (iDim = 0; iDim < nDim; iDim++) {
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            element->SetRef_Coord(iNode, iDim, val_Coord);
          }
        }

        /*--- In topology mode determine the penalty to apply to the mass,
         *    linear function of the physical density. ---*/
        su2double simp_penalty = 1.0;
        if (topology_mode) {
          simp_penalty = simp_minstiff+(1.0-simp_minstiff)*element_properties[iElem]->GetPhysicalDensity();
        }

        /*--- Set the properties of the element and compute its mass matrix. ---*/
        element->Set_ElProperties(element_properties[iElem]);

        numerics[FEA_TERM + thread*MAX_TERMS]->Compute_Mass_Matrix(element, config);

        /*--- Add contributions of this element to the mass matrix. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {
          for (jNode = 0; jNode < nNodes; jNode++) {

            su2double Mab = simp_penalty * element->Get_Mab(iNode, jNode);

            for (iVar = 0; iVar < nVar; iVar++) {
              TimeRes[indexNode[iNode]*nVar+iVar] += Mab * TimeRes_Aux.GetBlock(indexNode[iNode],iVar);
              TimeRes[indexNode[jNode]*nVar+iVar] += Mab * TimeRes_Aux.GetBlock(indexNode[jNode],iVar);
            }
          }
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  const bool prestretch_fem = config->GetPrestretch();

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {

          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

          for (iDim = 0; iDim < nDim; iDim++) {
            /*--- Compute current coordinate. ---*/
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            su2double val_Sol = nodes->GetSolution(indexNode[iNode],iDim) + val_Coord;

            /*--- If pre-stretched the reference coordinate is stored in the nodes. ---*/
            if (prestretch_fem)
              val_Coord = nodes->GetPrestretch(indexNode[iNode],iDim);

            /*--- Set coordinates. ---*/
            element->SetCurr_Coord(iNode, iDim, val_Sol);
            element->SetRef_Coord(iNode, iDim, val_Coord);
          }
        }

        /*--- In topology mode determine the penalty to apply to the stiffness ---*/
        su2double simp_penalty = 1.0;
        if (topology_mode) {
          su2double density = element_properties[iElem]->GetPhysicalDensity();
          simp_penalty = simp_minstiff+(1.0-simp_minstiff)*pow(density,simp_exponent);
        }

        /*--- Set the properties of the element. ---*/
        element->Set_ElProperties(element_properties[iElem]);

        /*--- Compute the components of the Jacobian and the stress term for the material. ---*/
        int NUM_TERM = thread*MAX_TERMS + element_based? element_properties[iElem]->GetMat_Mod() : FEA_TERM;

        numerics[NUM_TERM]->Compute_NodalStress_Term(element, config);

        for (iNode = 0; iNode < nNodes; iNode++) {
          auto Ta = element->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes[indexNode[iNode]*nVar+iVar] -= simp_penalty*Ta[iVar];
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_NodalStress(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  const bool prestretch_fem = config->GetPrestretch();

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();

  const unsigned short nStress = (nDim == 2) ? 3 : 6;

  /*--- Reduction variable to compute the maximum von Misses stress,
   *    each thread uses the first element of each row, the rest of
   *    the row is padding to prevent false sharing. ---*/
  su2activematrix auxMaxVonMisses(omp_get_max_threads(), 64/sizeof(su2double));
  auxMaxVonMisses.setConstant(0.0);

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Restart stress to avoid adding over results from previous time steps. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (unsigned short iStress = 0; iStress < nStress; iStress++) {
        nodes->SetStress_FEM(iPoint,iStress, 0.0);
      }
    }

#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, iDim, iVar, iStress;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {

          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

          for (iDim = 0; iDim < nDim; iDim++) {
            /*--- Compute current coordinate. ---*/
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            su2double val_Sol = nodes->GetSolution(indexNode[iNode],iDim) + val_Coord;

            /*--- If pre-stretched the reference coordinate is stored in the nodes. ---*/
            if (prestretch_fem)
              val_Coord = nodes->GetPrestretch(indexNode[iNode],iDim);

            /*--- Set coordinates. ---*/
            element->SetCurr_Coord(iNode, iDim, val_Sol);
            element->SetRef_Coord(iNode, iDim, val_Coord);
          }
        }

        /*--- In topology mode determine the penalty to apply to the stiffness ---*/
        su2double simp_penalty = 1.0;
        if (topology_mode) {
          simp_penalty = pow(element_properties[iElem]->GetPhysicalDensity(), simp_exponent);
        }

        /*--- Set the properties of the element. ---*/
        element->Set_ElProperties(element_properties[iElem]);

        /*--- Compute the averaged nodal stresses. ---*/
        int NUM_TERM = thread*MAX_TERMS + element_based? element_properties[iElem]->GetMat_Mod() : FEA_TERM;

        numerics[NUM_TERM]->Compute_Averaged_NodalStress(element, config);

        for (iNode = 0; iNode < nNodes; iNode++) {

          auto iPoint = indexNode[iNode];

          auto Ta = element->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysReact[iPoint*nVar+iVar] += simp_penalty*Ta[iVar];

          /*--- Divide the nodal stress by the number of elements that will contribute to this point. ---*/
          su2double weight = simp_penalty / geometry->node[iPoint]->GetnElem();

          for (iStress = 0; iStress < nStress; iStress++)
            nodes->AddStress_FEM(iPoint,iStress, weight*element->Get_NodalStress(iNode,iStress));
        }

      } // end iElem loop

    } // end color loop


    /*--- Compute the von Misses stress at each point, and the maximum for the domain. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {

      int thread = omp_get_thread_num();

      /*--- Get the stresses, added up from all the elements that connect to the node. ---*/

      auto Stress = nodes->GetStress_FEM(iPoint);
      su2double VonMises_Stress;

      if (nDim == 2) {

        su2double Sxx = Stress[0], Syy = Stress[1], Sxy = Stress[2];

        su2double S1, S2; S1 = S2 = (Sxx+Syy)/2;
        su2double tauMax = sqrt(pow((Sxx-Syy)/2, 2) + pow(Sxy,2));
        S1 += tauMax;
        S2 -= tauMax;

        VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);
      }
      else {

        su2double Sxx = Stress[0], Syy = Stress[1], Szz = Stress[3];
        su2double Sxy = Stress[2], Sxz = Stress[4], Syz = Stress[5];

        VonMises_Stress = sqrt(0.5*(pow(Sxx - Syy, 2) +
                                    pow(Syy - Szz, 2) +
                                    pow(Szz - Sxx, 2) +
                                    6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)));
      }

      nodes->SetVonMises_Stress(iPoint,VonMises_Stress);

      /*--- Update the maximum value of the Von Mises Stress ---*/

      auxMaxVonMisses(thread,0) = max(auxMaxVonMisses(thread,0), VonMises_Stress);
    }

  } // end SU2_OMP_PARALLEL

  /*--- Compute MaxVonMises_Stress across all threads and ranks. ---*/

  for(auto thread = 1ul; thread < auxMaxVonMisses.rows(); ++thread)
    auxMaxVonMisses(0,0) = max(auxMaxVonMisses(0,0), auxMaxVonMisses(thread,0));

  su2double MaxVonMises_Stress;
  SU2_MPI::Allreduce(auxMaxVonMisses.data(), &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/

  Total_CFEA = MaxVonMises_Stress;


  bool outputReactions = false;

  if (outputReactions) {

    bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

    ofstream myfile;
    myfile.open ("Reactions.txt");

    unsigned short iMarker, iDim, iVar;
    unsigned long iPoint, iVertex;
    su2double val_Reaction, val_Coord;

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
          for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
            for (iVar = 0; iVar < nVar; iVar++) {
              Residual[iVar] =   a_dt[0]*nodes->GetSolution_time_n(iPoint,iVar)        //a0*U(t)
              - a_dt[0]*nodes->GetSolution(iPoint,iVar)           //a0*U(t+dt)(k-1)
              + a_dt[2]*nodes->GetSolution_Vel_time_n(iPoint,iVar)    //a2*U'(t)
              + a_dt[3]*nodes->GetSolution_Accel_time_n(iPoint,iVar);  //a3*U''(t)
            }
            TimeRes_Aux.SetBlock(iPoint, Residual);
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
                    /*--- Retrieve the time contribution and reaction. ---*/
                    val_Reaction = LinSysReact.GetBlock(iPoint, iVar) + TimeRes.GetBlock(iPoint, iVar);
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

void CFEASolver::Compute_DeadLoad(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
#ifdef HAVE_OMP
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    auto chunkSize = roundUpDiv(OMP_MIN_SIZE, ColorGroupSize)*ColorGroupSize;

    /*--- Loop over element colors. ---*/
    for (auto color : ElemColoring)
    {
      SU2_OMP_FOR_DYN(chunkSize)
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];
#else
    /*--- Natural coloring. ---*/
    {
      for (auto iElem = 0ul; iElem < nElement; iElem++) {
#endif
        unsigned short iNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE];

        for (iNode = 0; iNode < nNodes; iNode++) {
          indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
          for (iDim = 0; iDim < nDim; iDim++) {
            su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
            element->SetRef_Coord(iNode, iDim, val_Coord);
          }
        }

        /*--- Penalize the dead load, do it by default to avoid unecessary "ifs", since it
         *    goes to the RHS there is no need to have a minimum value for stability ---*/
        su2double simp_penalty = element_properties[iElem]->GetPhysicalDensity();

        /*--- Set the properties of the element and compute its mass matrix. ---*/
        element->Set_ElProperties(element_properties[iElem]);

        numerics[FEA_TERM + thread*MAX_TERMS]->Compute_Dead_Load(element, config);

        /*--- Add contributions of this element to the mass matrix. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {

          auto Dead_Load = element->Get_FDL_a(iNode);

          su2double Aux_Dead_Load[MAXNVAR];
          for (iVar = 0; iVar < nVar; iVar++)
            Aux_Dead_Load[iVar] = simp_penalty*Dead_Load[iVar];

          nodes->Add_BodyForces_Res(indexNode[iNode], Aux_Dead_Load);
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

}

void CFEASolver::Compute_IntegrationConstants(CConfig *config) {

  su2double Delta_t= config->GetDelta_DynTime();

  su2double gamma = config->GetNewmark_gamma(), beta = config->GetNewmark_beta();

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      cout << "NOT IMPLEMENTED YET" << endl;
      break;
    case (NEWMARK_IMPLICIT):

      /*--- Integration constants for Newmark scheme ---*/

      a_dt[0]= 1 / (beta*pow(Delta_t,2.0));
      a_dt[1]= gamma / (beta*Delta_t);
      a_dt[2]= 1 / (beta*Delta_t);
      a_dt[3]= 1 /(2*beta) - 1;
      a_dt[4]= gamma/beta - 1;
      a_dt[5]= (Delta_t/2) * (gamma/beta - 2);
      a_dt[6]= Delta_t * (1-gamma);
      a_dt[7]= gamma * Delta_t;
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


void CFEASolver::BC_Clamped(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;

  bool dynamic = (config->GetTime_Domain());

  Solution[0] = 0.0;
  Solution[1] = 0.0;
  if (nDim==3) Solution[2] = 0.0;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Set and enforce solution at current and previous time-step ---*/
    nodes->SetSolution(iPoint,Solution);

    if (dynamic) {
      nodes->SetSolution_Vel(iPoint,Solution);
      nodes->SetSolution_Accel(iPoint,Solution);
      nodes->Set_Solution_time_n(iPoint,Solution);
      nodes->SetSolution_Vel_time_n(iPoint,Solution);
      nodes->SetSolution_Accel_time_n(iPoint,Solution);
    }

    /*--- Set and enforce 0 solution for mesh deformation ---*/
    nodes->SetBound_Disp(iPoint,Solution);

    LinSysSol.SetBlock(iPoint, Solution);
    LinSysReact.SetBlock(iPoint, Solution);
    Jacobian.EnforceSolutionAtNode(iPoint, Solution, LinSysRes);

  }

}

void CFEASolver::BC_Clamped_Post(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  bool dynamic = (config->GetTime_Domain());

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (nDim == 2) {
      Solution[0] = 0.0;  Solution[1] = 0.0;
    }
    else {
      Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
    }

    nodes->SetSolution(iPoint,Solution);

    if (dynamic) {
      nodes->SetSolution_Vel(iPoint,Solution);
      nodes->SetSolution_Accel(iPoint,Solution);
    }

  }

}

void CFEASolver::BC_DispDir(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iNode, iVertex;
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
  su2double ModAmpl = 1.0;

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

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Set and enforce solution ---*/
    LinSysSol.SetBlock(iNode, Disp_Dir);
    Jacobian.EnforceSolutionAtNode(iNode, Disp_Dir, LinSysRes);
  }

}

void CFEASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                                           unsigned short iMesh) {

  unsigned short iVar;
  unsigned long iPoint, total_index;

  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);    // Nonlinear analysis.
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);

  if (disc_adj_fem) {

    if (nonlinear_analysis) {

      /*--- For nonlinear discrete adjoint, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)): ABSOLUTE, norm of the incremental displacements ------------*/
      /*--- RTOL = norm(Residual(k): ABSOLUTE, norm of the residual (T-F) -----------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k): ABSOLUTE, norm of the product disp * res -----------*/

      Conv_Check[0] = LinSysSol.norm();               // Norm of the delta-solution vector
      Conv_Check[1] = LinSysRes.norm();               // Norm of the residual
      Conv_Check[2] = LinSysSol.dot(LinSysRes);       // Position for the energy tolerance

      /*--- MPI solution ---*/

      InitiateComms(geometry, config, SOLUTION_FEA);
      CompleteComms(geometry, config, SOLUTION_FEA);
    }
    else {
      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/
      /*---  Compute the residual Ax-f ---*/
#ifndef CODI_REVERSE_TYPE
      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
#else
      /*---  We need temporaries to interface with the matrix ---*/
      {
        CSysVector<passivedouble> sol, res;
        sol.PassiveCopy(LinSysSol);
        res.PassiveCopy(LinSysRes);
        CSysVector<passivedouble> aux(res);
        Jacobian.ComputeResidual(sol, res, aux);
        LinSysAux.PassiveCopy(aux);
      }
#endif

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
          AddRes_Max(iVar, fabs(LinSysAux[total_index]),
                     geometry->node[iPoint]->GetGlobalIndex(),
                     geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      InitiateComms(geometry, config, SOLUTION_FEA);
      CompleteComms(geometry, config, SOLUTION_FEA);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);

    }

  }
  else {
    if (nonlinear_analysis){

      /*--- MPI solution ---*/

      InitiateComms(geometry, config, SOLUTION_FEA);
      CompleteComms(geometry, config, SOLUTION_FEA);

    } else {

      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

      /*---  Compute the residual Ax-f ---*/
#ifndef CODI_REVERSE_TYPE
      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
#else
      /*---  We need temporaries to interface with the matrix ---*/
      {
        CSysVector<passivedouble> sol, res;
        sol.PassiveCopy(LinSysSol);
        res.PassiveCopy(LinSysRes);
        CSysVector<passivedouble> aux(res);
        Jacobian.ComputeResidual(sol, res, aux);
        LinSysAux.PassiveCopy(aux);
      }
#endif

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
          AddRes_Max(iVar, fabs(LinSysAux[total_index]),
                     geometry->node[iPoint]->GetGlobalIndex(),
                     geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      InitiateComms(geometry, config, SOLUTION_FEA);
      CompleteComms(geometry, config, SOLUTION_FEA);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);
    }

  }

}

void CFEASolver::BC_Normal_Displacement(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CFEASolver::BC_Normal_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Retrieve the normal pressure and the application conditions for the considered boundary ---*/

  su2double NormalLoad = config->GetLoad_Value(config->GetMarker_All_TagBound(val_marker));
  su2double TotalLoad = 0.0;

  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl = 1.0;

  su2double Ramp_Time = config->GetRamp_Time();

  ModAmpl = Compute_LoadCoefficient(CurrentTime, Ramp_Time, config);

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
          val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
          val_Sol = nodes->GetSolution(indexNode[iNode],iDim) + val_Coord;
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
        }
        else if (nonlinear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[1];
        }

        nodes->Add_SurfaceLoad_Res(indexNode[0],Residual);
        nodes->Add_SurfaceLoad_Res(indexNode[1],Residual);

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
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
          }

          nodes->Add_SurfaceLoad_Res(indexNode[0],Residual);
          nodes->Add_SurfaceLoad_Res(indexNode[1],Residual);
          nodes->Add_SurfaceLoad_Res(indexNode[2],Residual);

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
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
          }

          nodes->Add_SurfaceLoad_Res(indexNode[0],Residual);
          nodes->Add_SurfaceLoad_Res(indexNode[1],Residual);
          nodes->Add_SurfaceLoad_Res(indexNode[2],Residual);
          nodes->Add_SurfaceLoad_Res(indexNode[3],Residual);

        }

      }

    }

  }

}

void CFEASolver::BC_Dir_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

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
  su2double Ramp_Time = config->GetRamp_Time();

  su2double ModAmpl = 1.0;

  ModAmpl = Compute_LoadCoefficient(CurrentTime, Ramp_Time, config);

  TotalLoad = ModAmpl * LoadDirVal * LoadDirMult;

  /*--- Compute the norm of the vector that was passed in the config file ---*/
  su2double Norm = pow(Load_Dir_Local[0],2) + pow(Load_Dir_Local[1],2);
  if (nDim==3) Norm += pow(Load_Dir_Local[2],2);
  Norm = sqrt(Norm);
  
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

      nodes->Add_SurfaceLoad_Res(Point_0,Residual);
      nodes->Add_SurfaceLoad_Res(Point_1,Residual);

    }

    else {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {

        Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

        nodes->Add_SurfaceLoad_Res(Point_0,Residual);
        nodes->Add_SurfaceLoad_Res(Point_1,Residual);
        nodes->Add_SurfaceLoad_Res(Point_2,Residual);

      }
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {

        Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

        nodes->Add_SurfaceLoad_Res(Point_0,Residual);
        nodes->Add_SurfaceLoad_Res(Point_1,Residual);
        nodes->Add_SurfaceLoad_Res(Point_2,Residual);
        nodes->Add_SurfaceLoad_Res(Point_3,Residual);

      }

    }

  }

}

void CFEASolver::BC_Sine_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CFEASolver::BC_Damper(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

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

        dampValue = - 1.0 * dampC * nodes->GetSolution_Vel(Point_0,iVar);
        nodes->Set_SurfaceLoad_Res(Point_0,iVar, dampValue);

        dampValue = - 1.0 * dampC * nodes->GetSolution_Vel(Point_1,iVar);
        nodes->Set_SurfaceLoad_Res(Point_1,iVar, dampValue);
    }

    if (nDim == 3) {

      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);

      for (iVar = 0; iVar < nVar; iVar++){
          dampValue = - 1.0 * dampC * nodes->GetSolution_Vel(Point_2,iVar);
          nodes->Set_SurfaceLoad_Res(Point_2,iVar, dampValue);
      }

      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);
        for (iVar = 0; iVar < nVar; iVar++){
            dampValue = - 1.0 * dampC * nodes->GetSolution_Vel(Point_3,iVar);
            nodes->Set_SurfaceLoad_Res(Point_3,iVar, dampValue);
        }
      }

    }

  }

}

void CFEASolver::BC_Deforming(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker){

  unsigned short iDim;
  unsigned long iNode, iVertex;
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Retrieve the boundary displacement ---*/
    for (iDim = 0; iDim < nDim; iDim++) Solution[iDim] = nodes->GetBound_Disp(iNode,iDim);

    /*--- Set and enforce solution ---*/
    LinSysSol.SetBlock(iNode, Solution);
    Jacobian.EnforceSolutionAtNode(iNode, Solution, LinSysRes);

  }

}

void CFEASolver::Integrate_FSI_Loads(CGeometry *geometry, CConfig *config) {

  unsigned short iDim, iNode, nNode;
  unsigned long iPoint, iElem, nElem;

  unsigned short iMarkerInt, nMarkerInt = config->GetMarker_n_ZoneInterface()/2,
                 iMarker, nMarker = config->GetnMarker_All();

  /*--- Temporary storage to store the forces on the element faces ---*/
  vector<su2double> forces;

  /*--- Loop through the FSI interface pairs ---*/
  /*--- 1st pass to compute forces ---*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      /*--- Define the boundary element ---*/
      unsigned long nodeList[4];
      su2double coords[4][3];
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      for (iNode = 0; iNode < nNode; ++iNode) {
        nodeList[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; ++iDim)
          coords[iNode][iDim] = geometry->node[nodeList[iNode]]->GetCoord(iDim)+
                                nodes->GetSolution(nodeList[iNode],iDim);
      }

      /*--- Compute the area ---*/
      su2double area = 0.0;

      if (nDim == 2) {
        area = pow(coords[0][0]-coords[1][0],2) + pow(coords[0][1]-coords[1][1],2);
      }
      else {
        su2double a[3], b[3], Ni, Nj, Nk;

        if (!quad) { // sides of the triangle
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[1][iDim]-coords[0][iDim];
            b[iDim] = coords[2][iDim]-coords[0][iDim];
          }
        }
        else { // diagonals of the quadrilateral
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[2][iDim]-coords[0][iDim];
            b[iDim] = coords[3][iDim]-coords[1][iDim];
          }
        }
        /*--- Area = 0.5*||a x b|| ---*/
        Ni = a[1]*b[2]-a[2]*b[1];
        Nj =-a[0]*b[2]+a[2]*b[0];
        Nk = a[0]*b[1]-a[1]*b[0];

        area = 0.25*(Ni*Ni+Nj*Nj+Nk*Nk);
      }
      area = sqrt(area);

      /*--- Integrate ---*/
      passivedouble weight = 1.0/nNode;
      su2double force[3] = {0.0, 0.0, 0.0};

      for (iNode = 0; iNode < nNode; ++iNode)
        for (iDim = 0; iDim < nDim; ++iDim)
          force[iDim] += weight*area*nodes->Get_FlowTraction(nodeList[iNode],iDim);

      for (iDim = 0; iDim < nDim; ++iDim) forces.push_back(force[iDim]);
    }
  }

  /*--- 2nd pass to set values. This is to account for overlap in the markers. ---*/
  /*--- By putting the integrated values back into the nodes no changes have to be made elsewhere. ---*/
  nodes->Clear_FlowTraction();

  vector<su2double>::iterator force_it = forces.begin();

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;
      passivedouble weight = 1.0/nNode;

      su2double force[3];
      for (iDim = 0; iDim < nDim; ++iDim) force[iDim] = *(force_it++)*weight;

      for (iNode = 0; iNode < nNode; ++iNode) {
        iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
        nodes->Add_FlowTraction(iPoint,force);
      }
    }
  }

#ifdef HAVE_MPI
  /*--- Perform a global reduction, every rank will get the nodal values of all halo elements ---*/
  /*--- This should be cheaper than the "normal" way, since very few points are both halo and interface ---*/
  vector<unsigned long> halo_point_loc, halo_point_glb;
  vector<su2double> halo_force;

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      /*--- If this is an halo element we share the nodal forces ---*/
      for (iNode = 0; iNode < nNode; ++iNode)
        if (!geometry->node[geometry->bound[iMarker][iElem]->GetNode(iNode)]->GetDomain())
          break;

      if (iNode < nNode) {
        for (iNode = 0; iNode < nNode; ++iNode) {
          iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
          /*--- local is for when later we update the values in this rank ---*/
          halo_point_loc.push_back(iPoint);
          halo_point_glb.push_back(geometry->node[iPoint]->GetGlobalIndex());
          for (iDim = 0; iDim < nDim; ++iDim)
            halo_force.push_back(nodes->Get_FlowTraction(iPoint,iDim));
        }
      }
    }
  }
  /*--- Determine the size of the arrays we need ---*/
  unsigned long nHaloLoc = halo_point_loc.size();
  unsigned long nHaloMax;
  MPI_Allreduce(&nHaloLoc,&nHaloMax,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD);

  /*--- Shared arrays, all the: number of halo points; halo point global indices; respective forces ---*/
  unsigned long *halo_point_num = new unsigned long[size];
  unsigned long *halo_point_all = new unsigned long[size*nHaloMax];
  su2double *halo_force_all = new su2double[size*nHaloMax*nDim];

  /*--- Make "allgathers" extra safe by resizing all vectors to the same size (some
        issues observed when nHaloLoc = 0, especially with the discrete adjoint. ---*/
  halo_point_glb.resize(nHaloMax,0);
  halo_force.resize(nHaloMax*nDim,0.0);

  MPI_Allgather(&nHaloLoc,1,MPI_UNSIGNED_LONG,halo_point_num,1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  MPI_Allgather(halo_point_glb.data(),nHaloMax,MPI_UNSIGNED_LONG,halo_point_all,nHaloMax,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  SU2_MPI::Allgather(halo_force.data(),nHaloMax*nDim,MPI_DOUBLE,halo_force_all,nHaloMax*nDim,MPI_DOUBLE,MPI_COMM_WORLD);

  /*--- Find shared points with other ranks and update our values ---*/
  for (int proc = 0; proc < size; ++proc)
  if (proc != rank) {
    unsigned long offset = proc*nHaloMax;
    for (iPoint = 0; iPoint < halo_point_num[proc]; ++iPoint) {
      unsigned long iPoint_glb = halo_point_all[offset+iPoint];
      ptrdiff_t pos = find(halo_point_glb.begin(),halo_point_glb.end(),iPoint_glb)-halo_point_glb.begin();
      if (pos < long(halo_point_glb.size())) {
        unsigned long iPoint_loc = halo_point_loc[pos];
        nodes->Add_FlowTraction(iPoint_loc,&halo_force_all[(offset+iPoint)*nDim]);
      }
    }
  }

  delete [] halo_point_num;
  delete [] halo_point_all;
  delete [] halo_force_all;
#endif
}

su2double CFEASolver::Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, CConfig *config){

  su2double LoadCoeff = 1.0;

  bool Ramp_Load = config->GetRamp_Load();
  bool Sine_Load = config->GetSine_Load();
  bool Ramp_And_Release = config->GetRampAndRelease_Load();

  su2double SineAmp = 0.0, SineFreq = 0.0, SinePhase = 0.0;

  su2double TransferTime = 1.0;

  bool restart = config->GetRestart(); // Restart analysis
  bool fsi = config->GetFSI_Simulation();
  bool stat_fsi = !config->GetTime_Domain();

  /*--- This offset introduces the ramp load in dynamic cases starting from the restart point. ---*/
  bool offset = (restart && fsi && (!stat_fsi));
  su2double DeltaT = config->GetDelta_DynTime();
  su2double OffsetTime = 0.0;
  OffsetTime = DeltaT * (config->GetRestart_Iter()-1);

  /*--- Polynomial functions from https://en.wikipedia.org/wiki/Smoothstep ---*/

  if ((Ramp_Load) && (RampTime > 0.0)){

    TransferTime = CurrentTime / RampTime;

    if (offset) TransferTime = (CurrentTime - OffsetTime) / RampTime;

    switch (config->GetDynamic_LoadTransfer()) {
    case INSTANTANEOUS:
      LoadCoeff = 1.0;
      break;
    case POL_ORDER_1:
      LoadCoeff = TransferTime;
      break;
    case POL_ORDER_3:
      LoadCoeff = -2.0 * pow(TransferTime,3.0) + 3.0 * pow(TransferTime,2.0);
      break;
    case POL_ORDER_5:
      LoadCoeff = 6.0 * pow(TransferTime, 5.0) - 15.0 * pow(TransferTime, 4.0) + 10 * pow(TransferTime, 3.0);
      break;
    case SIGMOID_10:
      LoadCoeff = (1 / (1+exp(-1.0 * 10.0 * (TransferTime - 0.5)) ) );
      break;
    case SIGMOID_20:
      LoadCoeff = (1 / (1+exp(-1.0 * 20.0 * (TransferTime - 0.5)) ) );
      break;
    }

    if (TransferTime > 1.0) LoadCoeff = 1.0;

    LoadCoeff = max(LoadCoeff,0.0);
    LoadCoeff = min(LoadCoeff,1.0);

  }
  else if (Sine_Load){

      /*--- Retrieve amplitude, frequency (Hz) and phase (rad) ---*/
      SineAmp   = config->GetLoad_Sine()[0];
      SineFreq  = config->GetLoad_Sine()[1];
      SinePhase = config->GetLoad_Sine()[2];

      LoadCoeff = SineAmp * sin(2*PI_NUMBER*SineFreq*CurrentTime + SinePhase);
  }

  /*--- Add possibility to release the load after the ramp---*/
  if ((Ramp_And_Release) && (CurrentTime >  RampTime)){
    LoadCoeff = 0.0;
  }

  /*--- Store the force coefficient ---*/

  SetForceCoeff(LoadCoeff);

  return LoadCoeff;


}

void CFEASolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CFEASolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  const bool first_iter = (config->GetInnerIter() == 0);
  const bool dynamic = (config->GetTime_Domain());
  const bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);
  const bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);
  const bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);
  const bool body_forces = config->GetDeadLoad();

  /*--- For simplicity, no incremental loading is handled with increment of 1. ---*/
  const su2double loadIncr = config->GetIncrementalLoad()? loadIncrement : su2double(1.0);

  SU2_OMP_PARALLEL
  {
    unsigned long iPoint;
    unsigned short iVar;

    /*--- Loads common to static and dynamic problems. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- External surface load contribution. ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_SurfaceLoad_Res(iPoint,iVar);
      }

      /*--- Body forces contribution (dead load). ---*/

      if (body_forces) {
        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_BodyForces_Res(iPoint,iVar);
        }
      }

      /*--- FSI contribution (flow loads). ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_FlowTraction(iPoint,iVar);
      }

    }

    /*--- Dynamic contribution. ---*/

    if (dynamic) {

      /*--- Add the mass matrix contribution to the Jacobian. ---*/

      /*
       * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
       * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
       * iteration, as the Jacobian is recomputed.
       *
       * If the problem is linear, we add the Mass Matrix contribution to the Jacobian everytime because for
       * correct differentiation the Jacobian is recomputed every time step.
       *
       */
      if ((nonlinear_analysis && (newton_raphson || first_iter)) || linear_analysis) {
        Jacobian.MatrixMatrixAddition(a_dt[0], MassMatrix);
      }

      /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          TimeRes_Aux[iPoint*nVar+iVar] =
            a_dt[0]*nodes->GetSolution_time_n(iPoint,iVar) -      // a0*U(t)
            a_dt[0]*nodes->GetSolution(iPoint,iVar) +             // a0*U(t+dt)(k-1)
            a_dt[2]*nodes->GetSolution_Vel_time_n(iPoint,iVar) +  // a2*U'(t)
            a_dt[3]*nodes->GetSolution_Accel_time_n(iPoint,iVar); // a3*U''(t)
        }
      }

      /*--- Add M*TimeRes_Aux to the residual. ---*/

      MassMatrix.MatrixVectorProduct(TimeRes_Aux, TimeRes, geometry, config);
      LinSysRes += TimeRes;
    }

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  const bool dynamic = (config->GetTime_Domain());

  SU2_OMP_PARALLEL
  {
    unsigned long iPoint;
    unsigned short iVar;
  
    /*--- Update solution. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Displacements component of the solution. ---*/

        nodes->Add_DeltaSolution(iPoint, iVar, LinSysSol[iPoint*nVar+iVar]);

      }

    }

    if (dynamic) {
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        for (iVar = 0; iVar < nVar; iVar++) {

          /*--- Acceleration component of the solution ---*/
          /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

          su2double sol = a_dt[0]*(nodes->GetSolution(iPoint,iVar) -
                                   nodes->GetSolution_time_n(iPoint,iVar)) -
                          a_dt[2]* nodes->GetSolution_Vel_time_n(iPoint,iVar) -
                          a_dt[3]* nodes->GetSolution_Accel_time_n(iPoint,iVar);

          nodes->SetSolution_Accel(iPoint, iVar, sol);
        }

        for (iVar = 0; iVar < nVar; iVar++) {

          /*--- Velocity component of the solution ---*/
          /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

          su2double sol = nodes->GetSolution_Vel_time_n(iPoint,iVar)+
                          a_dt[6]* nodes->GetSolution_Accel_time_n(iPoint,iVar) +
                          a_dt[7]* nodes->GetSolution_Accel(iPoint,iVar);

          nodes->SetSolution_Vel(iPoint, iVar, sol);
        }
      }
    }

  } // end SU2_OMP_PARALLEL

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION_FEA);
  CompleteComms(geometry, config, SOLUTION_FEA);

}

void CFEASolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  const bool dynamic = (config->GetTime_Domain());

  SU2_OMP_PARALLEL
  {
    unsigned long iPoint;
    unsigned short iVar;

    /*--- Update solution and set it to be the solution after applying relaxation. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint=0; iPoint < nPointDomain; iPoint++) {
      nodes->SetSolution(iPoint, nodes->GetSolution_Pred(iPoint));
    }

    if (dynamic) {
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        for (iVar = 0; iVar < nVar; iVar++) {

          /*--- Acceleration component of the solution ---*/
          /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

          su2double sol = a_dt[0]*(nodes->GetSolution(iPoint,iVar) -
                                   nodes->GetSolution_time_n(iPoint,iVar)) -
                          a_dt[2]* nodes->GetSolution_Vel_time_n(iPoint,iVar) -
                          a_dt[3]* nodes->GetSolution_Accel_time_n(iPoint,iVar);

          nodes->SetSolution_Accel(iPoint, iVar, sol);
        }

        for (iVar = 0; iVar < nVar; iVar++) {

          /*--- Velocity component of the solution ---*/
          /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

          su2double sol = nodes->GetSolution_Vel_time_n(iPoint,iVar)+
                          a_dt[6]* nodes->GetSolution_Accel_time_n(iPoint,iVar) +
                          a_dt[7]* nodes->GetSolution_Accel(iPoint,iVar);

          nodes->SetSolution_Vel(iPoint, iVar, sol);
        }
      }
    }

    /*--- Perform the MPI communication of the solution ---*/
    SU2_OMP_MASTER
    {
      InitiateComms(geometry, config, SOLUTION_FEA);
      CompleteComms(geometry, config, SOLUTION_FEA);
    }
    SU2_OMP_BARRIER

    /*--- After the solution has been communicated, set the 'old' predicted solution as the solution. ---*/
    /*--- Loop over n points (as we have already communicated everything. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        nodes->SetSolution_Pred_Old(iPoint,iVar,nodes->GetSolution(iPoint,iVar));
      }
    }

  } // end SU2_OMP_PARALLEL

}


void CFEASolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  const bool first_iter = (config->GetInnerIter() == 0);
  const bool dynamic = (config->GetTime_Domain());
  const bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);
  const bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);
  const bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);
  const bool body_forces = config->GetDeadLoad();

  /*--- Blend between previous and current timestep. ---*/
  const su2double alpha_f = config->Get_Int_Coeffs(2);

  /*--- For simplicity, no incremental loading is handled with increment of 1. ---*/
  const su2double loadIncr = config->GetIncrementalLoad()? loadIncrement : su2double(1.0);

  SU2_OMP_PARALLEL
  {
    unsigned long iPoint;
    unsigned short iVar;

    /*--- Loads for static problems. ---*/
    if(!dynamic) {
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {

        /*--- External surface load contribution. ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_SurfaceLoad_Res(iPoint,iVar);
        }

        /*--- Body forces contribution (dead load). ---*/

        if (body_forces) {
          for (iVar = 0; iVar < nVar; iVar++) {
            LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_BodyForces_Res(iPoint,iVar);
          }
        }

        /*--- FSI contribution (flow loads). ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_FlowTraction(iPoint,iVar);
        }

      }
    }

    /*--- Loads for dynamic problems. ---*/

    if (dynamic) {

      /*--- Add the mass matrix contribution to the Jacobian. ---*/

      /*--- See notes on logic in ImplicitNewmark_Iteration(). ---*/
      if ((nonlinear_analysis && (newton_raphson || first_iter)) || linear_analysis) {
        Jacobian.MatrixMatrixAddition(a_dt[0], MassMatrix);
      }

      /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          TimeRes_Aux[iPoint*nVar+iVar] =
            a_dt[0]*nodes->GetSolution_time_n(iPoint,iVar) -      // a0*U(t)
            a_dt[0]*nodes->GetSolution(iPoint,iVar) +             // a0*U(t+dt)(k-1)
            a_dt[2]*nodes->GetSolution_Vel_time_n(iPoint,iVar) +  // a2*U'(t)
            a_dt[3]*nodes->GetSolution_Accel_time_n(iPoint,iVar); // a3*U''(t)
        }
      }

      /*--- Add M*TimeRes_Aux to the residual. ---*/

      MassMatrix.MatrixVectorProduct(TimeRes_Aux, TimeRes, geometry, config);
      LinSysRes += TimeRes;

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
        /*--- External surface load contribution ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes[iPoint*nVar+iVar] += loadIncr * ( (1-alpha_f) * nodes->Get_SurfaceLoad_Res(iPoint,iVar) +
                                                       alpha_f  * nodes->Get_SurfaceLoad_Res_n(iPoint,iVar) );
        }

        /*--- Add the contribution to the residual due to body forces.
         *--- It is constant over time, so it's not necessary to distribute it. ---*/

        if (body_forces) {
          for (iVar = 0; iVar < nVar; iVar++) {
            LinSysRes[iPoint*nVar+iVar] += loadIncr * nodes->Get_BodyForces_Res(iPoint,iVar);
          }
        }

        /*--- Add FSI contribution. ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes[iPoint*nVar+iVar] += loadIncr * ( (1-alpha_f) * nodes->Get_FlowTraction(iPoint,iVar) +
                                                       alpha_f  * nodes->Get_FlowTraction_n(iPoint,iVar) );
        }
      }
    }

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  /*--- Update displacement components of the solution. ---*/

  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      nodes->Add_DeltaSolution(iPoint, iVar, LinSysSol[iPoint*nVar+iVar]);

  /*--- Perform the MPI communication of the solution, displacements only. ---*/

  InitiateComms(geometry, config, SOLUTION_DISPONLY);
  CompleteComms(geometry, config, SOLUTION_DISPONLY);

}

void CFEASolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  const su2double alpha_f = config->Get_Int_Coeffs(2);
  const su2double alpha_m = config->Get_Int_Coeffs(3);

  /*--- Compute solution at t_n+1, and update velocities and accelerations ---*/

  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    unsigned short iVar;

    for (iVar = 0; iVar < nVar; iVar++) {

      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U(t+dt) = 1/alpha_f*(U(t+1-alpha_f)-alpha_f*U(t)) ---*/

      su2double sol = (1/(1-alpha_f)) * (nodes->GetSolution(iPoint,iVar) -
                                         alpha_f * nodes->GetSolution_time_n(iPoint,iVar));

      nodes->SetSolution(iPoint, iVar, sol);
    }

    for (iVar = 0; iVar < nVar; iVar++) {

      /*--- Acceleration component of the solution ---*/
      /*--- U''(t+dt-alpha_m) = a8*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

      su2double tmp = a_dt[8]*(nodes->GetSolution(iPoint,iVar) -
                               nodes->GetSolution_time_n(iPoint,iVar)) -
                      a_dt[2]* nodes->GetSolution_Vel_time_n(iPoint,iVar) -
                      a_dt[3]* nodes->GetSolution_Accel_time_n(iPoint,iVar);

      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U''(t+dt) = 1/alpha_m*(U''(t+1-alpha_m)-alpha_m*U''(t)) ---*/

      su2double sol = (1/(1-alpha_m)) * (tmp - alpha_m*nodes->GetSolution_Accel_time_n(iPoint,iVar));

      nodes->SetSolution_Accel(iPoint, iVar, sol);
    }

    for (iVar = 0; iVar < nVar; iVar++) {

      /*--- Velocity component of the solution ---*/
      /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

      su2double sol = nodes->GetSolution_Vel_time_n(iPoint,iVar)+
                      a_dt[6]* nodes->GetSolution_Accel_time_n(iPoint,iVar) +
                      a_dt[7]* nodes->GetSolution_Accel(iPoint,iVar);

      nodes->SetSolution_Vel(iPoint, iVar, sol);
    }

  }

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION_FEA);
  CompleteComms(geometry, config, SOLUTION_FEA);

}

void CFEASolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  /*--- Set the load conditions of the time step n+1 as the load conditions for time step n ---*/
  nodes->Set_SurfaceLoad_Res_n();
  nodes->Set_FlowTraction_n();

}

void CFEASolver::Solve_System(CGeometry *geometry, CConfig *config) {

  /*--- Initialize residual and solution at the ghost points ---*/

  for (auto iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

  /*--- Store the value of the residual. ---*/

  SetResLinSolver(System.GetResidual());
  
}


void CFEASolver::PredictStruct_Displacement(CGeometry **fea_geometry,
                                            CConfig *fea_config,
                                            CSolver ***fea_solution) {

  const unsigned short predOrder = fea_config->GetPredictorOrder();
  const su2double Delta_t = fea_config->GetDelta_DynTime();

  auto fea_nodes = fea_solution[MESH_0][FEA_SOL]->GetNodes();

  if(predOrder > 2 && rank == MASTER_NODE)
    cout << "Higher order predictor not implemented. Solving with order 0." << endl;

  /*--- To nPointDomain: we need to communicate the predicted solution after setting it. ---*/
  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint=0; iPoint < nPointDomain; iPoint++) {

    unsigned short iDim;

    switch (predOrder) {
      case 1: {
        const su2double* solDisp = fea_nodes->GetSolution(iPoint);
        const su2double* solVel = fea_nodes->GetSolution_Vel(iPoint);
        su2double* valPred = fea_nodes->GetSolution_Pred(iPoint);

        for (iDim=0; iDim < nDim; iDim++) {
          valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
        }
      } break;

      case 2: {
        const su2double* solDisp = fea_nodes->GetSolution(iPoint);
        const su2double* solVel = fea_nodes->GetSolution_Vel(iPoint);
        const su2double* solVel_tn = fea_nodes->GetSolution_Vel_time_n(iPoint);
        su2double* valPred = fea_nodes->GetSolution_Pred(iPoint);

        for (iDim=0; iDim < nDim; iDim++) {
          valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
        }
      } break;

      default: {
        fea_nodes->SetSolution_Pred(iPoint);
      } break;
    }

  }

}

void CFEASolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config,
                                           CSolver ***fea_solution, unsigned long iOuterIter) {

  unsigned long iPoint, iDim;
  su2double rbuf_numAitk = 0, sbuf_numAitk = 0;
  su2double rbuf_denAitk = 0, sbuf_denAitk = 0;

  const su2double *dispPred = nullptr;
  const su2double *dispCalc = nullptr;
  const su2double *dispPred_Old = nullptr;
  const su2double *dispCalc_Old = nullptr;
  su2double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
  su2double delta_deltaU[3] = {0.0, 0.0, 0.0};
  su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn_Min, WAitkDyn;

  unsigned short RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();

  /*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/

  if (RelaxMethod_FSI == NO_RELAXATION) {

    SetWAitken_Dyn(1.0);

  }
  else if (RelaxMethod_FSI == FIXED_PARAMETER) {

    SetWAitken_Dyn(fea_config->GetAitkenStatRelax());

  }
  else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {

    if (iOuterIter == 0) {

      WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
      WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();
      WAitkDyn_Min = fea_config->GetAitkenDynMinInit();

      WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);
      WAitkDyn = max(WAitkDyn, WAitkDyn_Min);

      SetWAitken_Dyn(WAitkDyn);

    }
    else {
      // To nPointDomain; we need to communicate the values
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        dispPred     = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution_Pred(iPoint);
        dispPred_Old = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution_Pred_Old(iPoint);
        dispCalc     = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution(iPoint);
        dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution_Old(iPoint);

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

      SU2_MPI::Allreduce(&sbuf_numAitk, &rbuf_numAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&sbuf_denAitk, &rbuf_denAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      WAitkDyn = GetWAitken_Dyn();

      if (rbuf_denAitk > EPS) {
        WAitkDyn = - 1.0 * WAitkDyn * rbuf_numAitk / rbuf_denAitk ;
      }

      WAitkDyn = max(WAitkDyn, 0.1);
      WAitkDyn = min(WAitkDyn, 1.0);

      SetWAitken_Dyn(WAitkDyn);

    }

  }
  else {
    if (rank == MASTER_NODE) cout << "No relaxation method used. " << endl;
  }

}

void CFEASolver::SetAitken_Relaxation(CGeometry **fea_geometry,
                                      CConfig *fea_config,
                                      CSolver ***fea_solution) {

  su2double WAitken = GetWAitken_Dyn();

  // To nPointDomain; we need to communicate the solutions (predicted, old and old predicted) after this routine
  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint=0; iPoint < nPointDomain; iPoint++) {

    /*--- Retrieve pointers to the predicted and calculated solutions ---*/
    su2double* dispPred = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution_Pred(iPoint);
    const su2double* dispCalc = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution(iPoint);

    /*--- Set predicted solution as the old predicted solution ---*/
    fea_solution[MESH_0][FEA_SOL]->GetNodes()->SetSolution_Pred_Old(iPoint);

    /*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
    fea_solution[MESH_0][FEA_SOL]->GetNodes()->SetSolution_Old(iPoint, dispCalc);

    /*--- Apply the Aitken relaxation ---*/
    for (unsigned short iDim=0; iDim < nDim; iDim++) {
      dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
    }
  }

}

void CFEASolver::Update_StructSolution(CGeometry **fea_geometry,
                                       CConfig *fea_config,
                                       CSolver ***fea_solution) {

  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint=0; iPoint < nPointDomain; iPoint++) {

    auto valSolutionPred = fea_solution[MESH_0][FEA_SOL]->GetNodes()->GetSolution_Pred(iPoint);

    fea_solution[MESH_0][FEA_SOL]->GetNodes()->SetSolution(iPoint, valSolutionPred);
  }

  /*--- Perform the MPI communication of the solution, displacements only ---*/

  InitiateComms(fea_geometry[MESH_0], fea_config, SOLUTION_DISPONLY);
  CompleteComms(fea_geometry[MESH_0], fea_config, SOLUTION_DISPONLY);

}


void CFEASolver::Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  unsigned long nTotalPoint = 1;

  bool dynamic = (config->GetTime_Domain());

  unsigned long TimeIter = config->GetTimeIter();

  su2double reference_geometry = 0.0, current_solution = 0.0;

  bool fsi = config->GetFSI_Simulation();

  su2double objective_function = 0.0, objective_function_reduce = 0.0;
  su2double weight_OF = 1.0;

  su2double objective_function_averaged = 0.0;

  SU2_MPI::Allreduce(&nPointDomain, &nTotalPoint, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  weight_OF = config->GetRefGeom_Penalty() / nTotalPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    for (iVar = 0; iVar < nVar; iVar++){

      /*--- Retrieve the value of the reference geometry ---*/
      reference_geometry = nodes->GetReference_Geometry(iPoint,iVar);

      /*--- Retrieve the value of the current solution ---*/
      current_solution = nodes->GetSolution(iPoint,iVar);

      /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
      objective_function += weight_OF * (current_solution - reference_geometry)*(current_solution - reference_geometry);
    }

  }

  SU2_MPI::Allreduce(&objective_function, &objective_function_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Total_OFRefGeom = objective_function_reduce + PenaltyValue;

  Global_OFRefGeom += Total_OFRefGeom;
  objective_function_averaged = Global_OFRefGeom / (TimeIter + 1.0 + EPS);

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
      averaged_gradient     = Total_ForwardGradient / (TimeIter + 1.0);
    }
    else {
      Total_ForwardGradient += local_forward_gradient;
      averaged_gradient      = Total_ForwardGradient / (TimeIter + 1.0);
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

void CFEASolver::Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;

  bool dynamic = config->GetTime_Domain();

  unsigned long TimeIter = config->GetTimeIter();

  su2double reference_geometry = 0.0, current_solution = 0.0;

  bool fsi = config->GetFSI_Simulation();

  su2double objective_function = 0.0, objective_function_reduce = 0.0;
  su2double distance_sq = 0.0 ;
  su2double weight_OF = 1.0;

  su2double objective_function_averaged = 0.0;

  /*--- TEMPORARY, for application in dynamic TestCase ---*/
  su2double difX = 0.0, difX_reduce = 0.0;
  su2double difY = 0.0, difY_reduce = 0.0;

  weight_OF = config->GetRefNode_Penalty();

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    if (geometry->node[iPoint]->GetGlobalIndex() == config->GetRefNode_ID() ){

      for (iVar = 0; iVar < nVar; iVar++){

        /*--- Retrieve the value of the reference geometry ---*/
        reference_geometry = config->GetRefNode_Displacement(iVar);

        /*--- Retrieve the value of the current solution ---*/
        current_solution = nodes->GetSolution(iPoint,iVar);

        /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
        distance_sq +=  (current_solution - reference_geometry)*(current_solution - reference_geometry);
      }

      objective_function = weight_OF * sqrt(distance_sq);

      difX = nodes->GetSolution(iPoint,0) - config->GetRefNode_Displacement(0);
      difY = nodes->GetSolution(iPoint,1) - config->GetRefNode_Displacement(1);

    }

  }

  SU2_MPI::Allreduce(&objective_function, &objective_function_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&difX, &difX_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&difY, &difY_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Total_OFRefNode = objective_function_reduce + PenaltyValue;

  Global_OFRefNode += Total_OFRefNode;
  objective_function_averaged = Global_OFRefNode / (TimeIter + 1.0 + EPS);

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
      averaged_gradient     = Total_ForwardGradient / (TimeIter + 1.0);
    }
    else {
      Total_ForwardGradient += local_forward_gradient;
      averaged_gradient      = Total_ForwardGradient / (TimeIter + 1.0);
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
      myfile_his << TimeIter << "\t";
      myfile_his << scientific << Total_OFRefNode << "\t";
      myfile_his << scientific << objective_function_averaged << "\t";
      myfile_his << scientific << difX_reduce << "\t";
      myfile_his << scientific << difY_reduce << endl;
      myfile_his.close();

    }

  }

}

void CFEASolver::Compute_OFVolFrac(CGeometry *geometry, CSolver **solver_container, CConfig *config)
{
  /*--- Perform a volume average of the physical density of the elements for topology optimization ---*/

  unsigned long iElem, nElem = geometry->GetnElem();
  su2double total_volume = 0.0, integral = 0.0, discreteness = 0.0;

  for (iElem=0; iElem<nElem; ++iElem) {
    /*--- count only elements that belong to the partition ---*/
    if ( geometry->node[geometry->elem[iElem]->GetNode(0)]->GetDomain() ){
      su2double volume = geometry->elem[iElem]->GetVolume();
      su2double rho = element_properties[iElem]->GetPhysicalDensity();
      total_volume += volume;
      integral += volume*rho;
      discreteness += volume*4.0*rho*(1.0-rho);
    }
  }

  su2double tmp;
  SU2_MPI::Allreduce(&total_volume,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  total_volume = tmp;
  SU2_MPI::Allreduce(&integral,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  integral = tmp;
  SU2_MPI::Allreduce(&discreteness,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  discreteness = tmp;

  if (config->GetKind_ObjFunc() == TOPOL_DISCRETENESS)
    Total_OFVolFrac = discreteness/total_volume;
  else
    Total_OFVolFrac = integral/total_volume;

  // TODO: Temporary output file for the objective function. Will be integrated in the output once is refurbished.
  if (rank == MASTER_NODE){
    cout << "Objective function: " << Total_OFVolFrac << "." << endl;

    ofstream myfile_res;
    if (config->GetKind_ObjFunc() == TOPOL_DISCRETENESS)
      myfile_res.open ("of_topdisc.dat");
    else
      myfile_res.open ("of_volfrac.dat");

    myfile_res.precision(15);
    myfile_res << scientific << Total_OFVolFrac << endl;
    myfile_res.close();
  }
}

void CFEASolver::Compute_OFCompliance(CGeometry *geometry, CSolver **solver_container, CConfig *config)
{
  unsigned long iPoint;
  unsigned short iVar;
  su2double nodalForce[MAXNVAR];

  /*--- Types of loads to consider ---*/
  bool body_forces = config->GetDeadLoad();

  /*--- If the loads are being applied incrementaly ---*/
  bool incremental_load = config->GetIncrementalLoad();

  /*--- Computation (compliance = sum(f dot u) ) ---*/
  /*--- Cannot be computed as u^T K u as the problem may not be linear ---*/

  Total_OFCompliance = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Initialize with loads speficied through config ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      nodalForce[iVar] = nodes->Get_SurfaceLoad_Res(iPoint,iVar);

    /*--- Add contributions due to body forces ---*/
    if (body_forces)
      for (iVar = 0; iVar < nVar; iVar++)
        nodalForce[iVar] += nodes->Get_BodyForces_Res(iPoint,iVar);

    /*--- Add contributions due to fluid loads---*/
    for (iVar = 0; iVar < nVar; iVar++)
      nodalForce[iVar] += nodes->Get_FlowTraction(iPoint,iVar);

    /*--- Correct for incremental loading ---*/
    if (incremental_load)
      for (iVar = 0; iVar < nVar; iVar++)
        nodalForce[iVar] *= loadIncrement;

    /*--- Add work contribution from this node ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Total_OFCompliance += nodalForce[iVar]*nodes->GetSolution(iPoint,iVar);
  }

  su2double tmp;
  SU2_MPI::Allreduce(&Total_OFCompliance,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  Total_OFCompliance = tmp;

  // TODO: Temporary output file for the objective function. Will be integrated in the output once refurbished.
  if (rank == MASTER_NODE) {
    cout << "Objective function: " << Total_OFCompliance << "." << endl;

    ofstream file;
    file.open("of_topcomp.dat");
    file.precision(15);
    file << scientific << Total_OFCompliance << endl;
    file.close();
  }
}

void CFEASolver::Stiffness_Penalty(CGeometry *geometry, CSolver **solver, CNumerics **numerics, CConfig *config){

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

    switch(geometry->elem[iElem]->GetVTK_Type()) {
      case TRIANGLE:      nNodes = 3; EL_KIND = EL_TRIA; break;
      case QUADRILATERAL: nNodes = 4; EL_KIND = EL_QUAD; break;
      case TETRAHEDRON:   nNodes = 4; EL_KIND = EL_TETRA; break;
      case PYRAMID:       nNodes = 5; EL_KIND = EL_PYRAM; break;
      case PRISM:         nNodes = 6; EL_KIND = EL_PRISM; break;
      case HEXAHEDRON:    nNodes = 8; EL_KIND = EL_HEXA; break;
    }

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; iDim++) {
          val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
          val_Sol = nodes->GetSolution(indexNode[iNode],iDim) + val_Coord;
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
          element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(iNode, iDim, val_Sol);
        }
    }

    // Avoid double-counting elements:
    // Only add the value if the first node is in the domain
    if (geometry->node[indexNode[0]]->GetDomain()) {

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

  SU2_MPI::Allreduce(&weightedValue, &weightedValue_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&totalVolume, &totalVolume_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  ratio = 1.0 - weightedValue_reduce/totalVolume_reduce;

  PenaltyValue = config->GetTotalDV_Penalty() * ratio * ratio;

}

void CFEASolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  unsigned short iVar, nSolVar;
  unsigned long index;

  string filename;

  bool dynamic = (config->GetTime_Domain());
  bool fluid_structure = config->GetFSI_Simulation();
  bool discrete_adjoint = config->GetDiscrete_Adjoint();

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

  su2double *Sol = new su2double[nSolVar];

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Restart the solution from file information ---*/

  filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

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
        nodes->SetSolution(iPoint_Local,iVar, Sol[iVar]);
        if (dynamic) {
          nodes->Set_Solution_time_n(iPoint_Local,iVar, Sol[iVar]);
          nodes->SetSolution_Vel(iPoint_Local,iVar, Sol[iVar+nVar]);
          nodes->SetSolution_Vel_time_n(iPoint_Local,iVar, Sol[iVar+nVar]);
          nodes->SetSolution_Accel(iPoint_Local,iVar, Sol[iVar+2*nVar]);
          nodes->SetSolution_Accel_time_n(iPoint_Local,iVar, Sol[iVar+2*nVar]);
        }
        if (fluid_structure && !dynamic) {
          nodes->SetSolution_Pred(iPoint_Local,iVar, Sol[iVar]);
          nodes->SetSolution_Pred_Old(iPoint_Local,iVar, Sol[iVar]);
        }
        if (fluid_structure && discrete_adjoint){
          nodes->SetSolution_Old(iPoint_Local,iVar, Sol[iVar]);
        }
      }
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI. If dynamic, we also need to communicate the old solution ---*/

  solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_FEA);
  solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_FEA);

  if (dynamic) {
    solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_FEA_OLD);
    solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_FEA_OLD);
  }
  if (fluid_structure && !dynamic){
      solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_PRED);
      solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_PRED);

      solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_PRED_OLD);
      solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_PRED_OLD);
  }

  delete [] Sol;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CFEASolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset)
{
  /*--- Register the element density to get the derivatives required for
  material-based topology optimization, this is done here because element_properties
  is a member of CFEASolver only. ---*/
  if (!config->GetTopology_Optimization()) return;

  for (unsigned long iElem = 0; iElem < geometry->GetnElem(); iElem++)
    element_properties[iElem]->RegisterDensity();
}

void CFEASolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config)
{
  /*--- Extract and output derivatives for topology optimization, this is done
  here because element_properties is a member of CFEASolver only and the output
  structure only supports nodal values (these are elemental). ---*/
  if (!config->GetTopology_Optimization()) return;

  unsigned long iElem,
                nElem = geometry->GetnElem(),
                nElemDomain = geometry->GetGlobal_nElemDomain();

  /*--- Allocate and initialize an array onto which the derivatives of every partition
  will be reduced, this is to output results in the correct order, it is not a very
  memory efficient solution... single precision is enough for output. ---*/
  float *send_buf = new float[nElemDomain], *rec_buf = nullptr;
  for(iElem=0; iElem<nElemDomain; ++iElem) send_buf[iElem] = 0.0;

  for(iElem=0; iElem<nElem; ++iElem) {
    unsigned long iElem_global = geometry->elem[iElem]->GetGlobalIndex();
    send_buf[iElem_global] = SU2_TYPE::GetValue(element_properties[iElem]->GetAdjointDensity());
  }

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) rec_buf = new float[nElemDomain];
  /*--- Need to use this version of Reduce instead of the wrapped one because we use float ---*/
  MPI_Reduce(send_buf,rec_buf,nElemDomain,MPI_FLOAT,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
#else
  rec_buf = send_buf;
#endif

  /*--- The master writes the file ---*/
  if (rank == MASTER_NODE) {
    string filename = config->GetTopology_Optim_FileName();
    ofstream file;
    file.open(filename.c_str());
    for(iElem=0; iElem<nElemDomain; ++iElem) file << rec_buf[iElem] << "\n";
    file.close();
  }

  delete [] send_buf;
#ifdef HAVE_MPI
  if (rank == MASTER_NODE) delete [] rec_buf;
#endif

}

void CFEASolver::FilterElementDensities(CGeometry *geometry, CConfig *config)
{
  /*--- Apply a filter to the design densities of the elements to generate the
  physical densities which are the ones used to penalize their stiffness. ---*/

  unsigned short type, search_lim;
  su2double param, radius;

  vector<pair<unsigned short,su2double> > kernels;
  vector<su2double> filter_radius;
  for (unsigned short iKernel=0; iKernel<config->GetTopology_Optim_Num_Kernels(); ++iKernel)
  {
    config->GetTopology_Optim_Kernel(iKernel,type,param,radius);
    kernels.push_back(make_pair(type,param));
    filter_radius.push_back(radius);
  }
  search_lim = config->GetTopology_Search_Limit();

  unsigned long iElem, nElem = geometry->GetnElem();

  su2double *design_rho = new su2double [nElem],
            *physical_rho = new su2double [nElem];

  /*--- "Rectify" the input ---*/
  for (iElem=0; iElem<nElem; ++iElem) {
    su2double rho = element_properties[iElem]->GetDesignDensity();
    if      (rho > 1.0) design_rho[iElem] = 1.0;
    else if (rho < 0.0) design_rho[iElem] = 0.0;
    else                design_rho[iElem] = rho;
  }

  geometry->FilterValuesAtElementCG(filter_radius, kernels, search_lim, design_rho, physical_rho);

  /*--- Apply projection ---*/
  config->GetTopology_Optim_Projection(type,param);
  switch (type) {
    case NO_PROJECTION: break;
    case HEAVISIDE_UP:
      for (iElem=0; iElem<nElem; ++iElem)
        physical_rho[iElem] = 1.0-exp(-param*physical_rho[iElem])+physical_rho[iElem]*exp(-param);
      break;
    case HEAVISIDE_DOWN:
      for (iElem=0; iElem<nElem; ++iElem)
        physical_rho[iElem] = exp(-param*(1.0-physical_rho[iElem]))-(1.0-physical_rho[iElem])*exp(-param);
      break;
    default:
      SU2_MPI::Error("Unknown type of projection function",CURRENT_FUNCTION);
  }

  /*--- If input was out of bounds use the bound instead of the filtered
   value, useful to enforce solid or void regions (e.g. a skin). ---*/
  for (iElem=0; iElem<nElem; ++iElem) {
    su2double rho = element_properties[iElem]->GetDesignDensity();
    if      (rho > 1.0) element_properties[iElem]->SetPhysicalDensity(1.0);
    else if (rho < 0.0) element_properties[iElem]->SetPhysicalDensity(0.0);
    else element_properties[iElem]->SetPhysicalDensity(physical_rho[iElem]);
  }

  delete [] design_rho;
  delete [] physical_rho;

  /*--- For when this method is called directly, e.g. by the adjoint solver. ---*/
  topol_filter_applied = true;
}
