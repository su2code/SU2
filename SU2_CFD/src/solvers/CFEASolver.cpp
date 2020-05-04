/*!
 * \file CFEASolver.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
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

#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/variables/CFEABoundVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

using namespace GeometryToolbox;


CFEASolver::CFEASolver(bool mesh_deform_mode) : CSolver(mesh_deform_mode) {

  nElement = 0;
  nDim = 0;
  nMarker = 0;

  nPoint = 0;
  nPointDomain = 0;

  Total_CFEA = 0.0;
  WAitken_Dyn = 0.0;
  WAitken_Dyn_tn1 = 0.0;

  element_container = new CElement** [MAX_TERMS]();
  for (unsigned short iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS*omp_get_max_threads()]();

  topol_filter_applied = false;
  element_based = false;
  initial_calc = true;

}

CFEASolver::CFEASolver(CGeometry *geometry, CConfig *config) : CSolver() {

  unsigned short iVar;

  bool dynamic = (config->GetTime_Domain());

  /*--- Test whether we consider dielectric elastomers ---*/
  bool de_effects = config->GetDE_Effects();

  /*--- A priori we don't have an element-based input file (most of the applications will be like this) ---*/
  element_based = false;
  topol_filter_applied = false;
  initial_calc = true;

  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();

  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();

  /*--- Here is where we assign the kind of each element ---*/

  /*--- First level: different possible terms of the equations ---*/
  element_container = new CElement** [MAX_TERMS]();
  for (unsigned short iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS*omp_get_max_threads()]();

  SU2_OMP_PARALLEL
  {
    const int offset = omp_get_thread_num()*MAX_FE_KINDS;

    if (nDim == 2) {
      /*--- Basic terms ---*/
      element_container[FEA_TERM][EL_TRIA+offset] = new CTRIA1();
      element_container[FEA_TERM][EL_QUAD+offset] = new CQUAD4();

      if (de_effects) {
        element_container[DE_TERM][EL_TRIA+offset] = new CTRIA1();
        element_container[DE_TERM][EL_QUAD+offset] = new CQUAD4();
      }
    }
    else {
      element_container[FEA_TERM][EL_TETRA+offset] = new CTETRA1();
      element_container[FEA_TERM][EL_HEXA +offset] = new CHEXA8 ();
      element_container[FEA_TERM][EL_PYRAM+offset] = new CPYRAM5();
      element_container[FEA_TERM][EL_PRISM+offset] = new CPRISM6();

      if (de_effects) {
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

  SetFSI_ConvValue(0,0.0);
  SetFSI_ConvValue(1,0.0);

  nVar = nDim;

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS = new su2double[nVar]();
  Residual_Max = new su2double[nVar]();
  Point_Max = new unsigned long[nVar]();
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim]();
  }

  /*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/

  unsigned short nSolVar;
  string text_line, filename;
  ifstream restart_file;

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

  su2double* SolInit = new su2double[nSolVar]();

  /*--- Initialize from zero everywhere ---*/

  nodes = new CFEABoundVariable(SolInit, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  delete [] SolInit;

  /*--- Set which points are vertices and allocate boundary data. ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
      if (iVertex >= 0) {
        nodes->Set_isVertex(iPoint,true);
        break;
      }
    }
  static_cast<CFEABoundVariable*>(nodes)->AllocateBoundaryVariables(config);


  if (config->GetRefGeom()) Set_ReferenceGeometry(geometry, config);
  if (config->GetPrestretch()) Set_Prestretch(geometry, config);

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
  LinSysReact.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Initialize structures for hybrid-parallel mode. ---*/
  HybridParallelInitialization(geometry);

  /*--- Initialize the value of the total objective function ---*/
  Total_OFRefGeom = 0.0;
  Total_OFRefNode = 0.0;
  Total_OFVolFrac = 0.0;

  /*--- Initialize the value of the global objective function ---*/
  Global_OFRefGeom = 0.0;
  Global_OFRefNode = 0.0;

  /*--- Initialize the value of the total gradient for the forward mode ---*/
  Total_ForwardGradient = 0.0;

  OutputForwardModeGradient(config, true, 0.0, 0.0, 0.0, 0.0);

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
  else {
    ForceCoeff = 1.0;
  }

  /*--- Penalty value - to maintain constant the stiffness in optimization problems - TODO: this has to be improved ---*/
  PenaltyValue = 0.0;

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION_FEA);
  CompleteComms(geometry, config, SOLUTION_FEA);

  /*--- If dynamic, we also need to communicate the old solution ---*/

  if (dynamic) {
    InitiateComms(geometry, config, SOLUTION_FEA_OLD);
    CompleteComms(geometry, config, SOLUTION_FEA_OLD);
  }

  if (size != SINGLE_NODE) {
    vector<unsigned short> essentialMarkers;
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      const auto kindBnd = config->GetMarker_All_KindBC(iMarker);
      if ((kindBnd == CLAMPED_BOUNDARY) ||
          (kindBnd == DISP_DIR_BOUNDARY) ||
          (kindBnd == DISPLACEMENT_BOUNDARY)) {
        essentialMarkers.push_back(iMarker);
      }
    }
    Set_VertexEliminationSchedule(geometry, essentialMarkers);
  }

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "FEA";
}

CFEASolver::~CFEASolver(void) {

  if (element_container != nullptr) {
    for (unsigned int iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (unsigned int jVar = 0; jVar < MAX_FE_KINDS*omp_get_max_threads(); jVar++) {
        delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }

  if (element_properties != nullptr) {
    for (unsigned long iElem = 0; iElem < nElement; iElem++)
      delete element_properties[iElem];
    delete [] element_properties;
  }

  delete [] iElem_iDe;

  delete nodes;

  if (LockStrategy) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      omp_destroy_lock(&UpdateLocks[iPoint]);
  }
}

void CFEASolver::HybridParallelInitialization(CGeometry* geometry) {
#ifdef HAVE_OMP
  /*--- Get the element coloring. ---*/

  su2double parallelEff = 1.0;
  const auto& coloring = geometry->GetElementColoring(&parallelEff);

  /*--- If the coloring is too bad use lock-guarded accesses
   *    to CSysMatrix/Vector in element loops instead. ---*/
  LockStrategy = parallelEff < COLORING_EFF_THRESH;

  /*--- When using locks force a single color to reduce the color loop overhead. ---*/
  if (LockStrategy && (coloring.getOuterSize()>1))
    geometry->SetNaturalElementColoring();

  if (!coloring.empty()) {
    /*--- We are not constrained by the color group size when using locks. ---*/
    auto groupSize = LockStrategy? 1ul : geometry->GetElementColorGroupSize();
    auto nColor = coloring.getOuterSize();
    ElemColoring.reserve(nColor);

    for(auto iColor = 0ul; iColor < nColor; ++iColor)
      ElemColoring.emplace_back(coloring.innerIdx(iColor), coloring.getNumNonZeros(iColor), groupSize);
  }

  su2double minEff = 1.0;
  SU2_MPI::Reduce(&parallelEff, &minEff, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);

  if (minEff < COLORING_EFF_THRESH) {
    cout << "WARNING: The element coloring efficiency was " << minEff << ", a fallback strategy is in use.\n"
         << "         Better performance may be possible by reducing the number of threads per rank." << endl;
  }

  if (LockStrategy) {
    UpdateLocks.resize(nPoint);
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      omp_init_lock(&UpdateLocks[iPoint]);
  }

  omp_chunk_size = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_SIZE);
#else
  ElemColoring[0] = DummyGridColor<>(nElement);
#endif
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
      element_properties[iElem] = new CElementProperty(FEA_TERM, 0, 0, 0);
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

        if (config->GetAdvanced_FEAElementBased() || topology_mode){
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];
          element_properties[iElem_Local] = new CElementProperty(
            elProperties[0], elProperties[1], elProperties[2], elProperties[3]);
        }
        else{
          point_line >> index >> elProperties[0];
          element_properties[iElem_Local] = new CElementProperty(0, elProperties[0], 0, 0);
        }


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

    if (iElem_Global_Local != nElement) {
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

  /*--- The first line is the header ---*/

  getline (prestretch_file, text_line);

  while (getline (prestretch_file, text_line)) {
    istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    MI = Global2Local.find(iPoint_Global);
    if (MI != Global2Local.end()) {

      iPoint_Local = Global2Local[iPoint_Global];

      su2double Sol[MAXNVAR] = {0.0};

      if (nDim == 2) point_line >> Sol[0] >> Sol[1] >> index;
      if (nDim == 3) point_line >> Sol[0] >> Sol[1] >> Sol[2] >> index;

      for (iVar = 0; iVar < nVar; iVar++) nodes->SetPrestretch(iPoint_Local,iVar, Sol[iVar]);

      iPoint_Global_Local++;
    }
    iPoint_Global++;
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local != nPointDomain) {
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
  su2double *Buffer_Receive_U = nullptr, *Buffer_Send_U = nullptr;

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
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, nullptr);

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          nodes->SetPrestretch(iPoint, iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

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

      su2double Sol[MAXNVAR] = {0.0};

      if (nDim == 2){
        Sol[0] = PrintingToolbox::stod(point_line[3]);
        Sol[1] = PrintingToolbox::stod(point_line[4]);
      } else {
        Sol[0] = PrintingToolbox::stod(point_line[4]);
        Sol[1] = PrintingToolbox::stod(point_line[5]);
        Sol[2] = PrintingToolbox::stod(point_line[6]);
      }

      for (iVar = 0; iVar < nVar; iVar++)
        nodes->SetReference_Geometry(iPoint_Local, iVar, Sol[iVar]);

      iPoint_Global_Local++;
    }
    iPoint_Global++;
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n")  +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Close the restart file ---*/

  reference_file.close();

  /*--- Free memory needed for the transformation ---*/

  delete [] Global2Local;

}

void CFEASolver::Set_VertexEliminationSchedule(CGeometry *geometry, const vector<unsigned short>& markers) {

  /*--- Store global point indices of essential BC markers. ---*/
  vector<unsigned long> myPoints;

  for (auto iMarker : markers) {
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      myPoints.push_back(geometry->node[iPoint]->GetGlobalIndex());
    }
  }

  const unordered_set<unsigned long> markerPoints(myPoints.begin(), myPoints.end());

  vector<unsigned long> numPoints(size);
  unsigned long num = myPoints.size();
  SU2_MPI::Allgather(&num, 1, MPI_UNSIGNED_LONG, numPoints.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Global to local map for the halo points of the rank (not covered by the CGeometry map). ---*/
  unordered_map<unsigned long, unsigned long> Global2Local;
  for (auto iPoint = nPointDomain; iPoint < nPoint; ++iPoint) {
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
  }

  /*--- Populate elimination list. ---*/
  ExtraVerticesToEliminate.clear();

  for (int i = 0; i < size; ++i) {
    /*--- Send our point list. ---*/
    if (rank == i) {
      SU2_MPI::Bcast(myPoints.data(), numPoints[i], MPI_UNSIGNED_LONG, rank, MPI_COMM_WORLD);
      continue;
    }

    /*--- Receive point list. ---*/
    vector<unsigned long> theirPoints(numPoints[i]);
    SU2_MPI::Bcast(theirPoints.data(), numPoints[i], MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);

    for (auto iPointGlobal : theirPoints) {
      /*--- Check if the rank has the point. ---*/
      auto it = Global2Local.find(iPointGlobal);
      if (it == Global2Local.end()) continue;

      /*--- If the point is not covered by this rank's markers, mark it for elimination. ---*/
      if (markerPoints.count(iPointGlobal) == 0)
        ExtraVerticesToEliminate.push_back(it->second);
    }
  }

}

void CFEASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics,
                               unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {

  const bool dynamic = config->GetTime_Domain();
  const bool first_iter = (config->GetInnerIter() == 0);
  const bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  const bool body_forces = config->GetDeadLoad();
  const bool fsi = config->GetFSI_Simulation();
  const bool topology_mode = config->GetTopology_Optimization();

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
    FilterElementDensities(geometry, config);
    topol_filter_applied = true;
  }

  /*
   * If the problem is dynamic we need a mass matrix, which will be constant along the calculation
   * both for linear and nonlinear analysis, i.e. only computed once on the first time step.
   *
   * The same with the integration constants, as for now we consider the time step to be constant.
   */
  if (dynamic && (initial_calc || disc_adj_fem)) {
    Compute_IntegrationConstants(config);
    Compute_MassMatrix(geometry, numerics, config);
  }

  /*
   * If body forces are taken into account, we need to compute the term that goes into the residual,
   * which will be constant along the calculation both for linear and nonlinear analysis.
   *
   * Only initialized once, at the first iteration of the first time step.
   */
  if (body_forces && (initial_calc || disc_adj_fem))
    Compute_DeadLoad(geometry, numerics, config);

  /*--- Clear the linear system solution. ---*/
  SU2_OMP_PARALLEL
  {
    LinSysSol.SetValZero();
  }

  /*--- Clear external forces. ---*/
  nodes->Clear_SurfaceLoad_Res();

  /*
   * FSI loads (computed upstream) need to be integrated if a nonconservative interpolation scheme is in use
   */
  if (fsi && first_iter && (idxIncrement==0)) Integrate_FSI_Loads(geometry, config);

  /*--- Next call to Preprocessing will not be "initial_calc" and linear operations will not be repeated. ---*/
  initial_calc = false;

}

void CFEASolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  SU2_OMP_PARALLEL
  {
  if (!config->GetPrestretch()) {
    su2double zeros[MAXNVAR] = {0.0};
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
      nodes->SetSolution(iPoint, zeros);
  }
  else {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
      nodes->SetSolution(iPoint, nodes->GetPrestretch(iPoint));
  }
  } // end parallel
}

void CFEASolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Clear vector and matrix before calculation. ---*/
    LinSysRes.SetValZero();
    Jacobian.SetValZero();

    for(auto color : ElemColoring) {

      /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
      SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];

        unsigned short iNode, jNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE_3D];

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
        int NUM_TERM = thread*MAX_TERMS + element_properties[iElem]->GetMat_Mod();

        numerics[NUM_TERM]->Compute_Tangent_Matrix(element, config);

        /*--- Update residual and stiffness matrix with contributions from the element. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {

          if (LockStrategy) omp_set_lock(&UpdateLocks[indexNode[iNode]]);

          auto Ta = element->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes(indexNode[iNode], iVar) -= simp_penalty*Ta[iVar];

          for (jNode = 0; jNode < nNodes; jNode++) {
            auto Kab = element->Get_Kab(iNode, jNode);
            Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], simp_penalty, Kab);
          }

          if (LockStrategy) omp_unset_lock(&UpdateLocks[indexNode[iNode]]);
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  const bool prestretch_fem = config->GetPrestretch();
  const bool de_effects = config->GetDE_Effects();

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Clear vector and matrix before calculation. ---*/
    LinSysRes.SetValZero();
    Jacobian.SetValZero();

    for(auto color : ElemColoring) {

      /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
      SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];

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
        unsigned long indexNode[MAXNNODE_3D];

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
        int NUM_TERM = thread*MAX_TERMS + element_properties[iElem]->GetMat_Mod();

        numerics[NUM_TERM]->Compute_Tangent_Matrix(fea_elem, config);

        /*--- Compute the electric component of the Jacobian and the stress term. ---*/
        if (de_effects)
          numerics[DE_TERM + thread*MAX_TERMS]->Compute_Tangent_Matrix(de_elem, config);

        /*--- Update residual and stiffness matrix with contributions from the element. ---*/
        for (iNode = 0; iNode < nNodes; iNode++) {

          if (LockStrategy) omp_set_lock(&UpdateLocks[indexNode[iNode]]);

          auto Ta = fea_elem->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes(indexNode[iNode], iVar) -= simp_penalty*Ta[iVar];

          /*--- Retrieve the electric contribution to the residual. ---*/
          if (de_effects) {
            auto Ta_DE = de_elem->Get_Kt_a(iNode);
            for (iVar = 0; iVar < nVar; iVar++)
              LinSysRes(indexNode[iNode], iVar) -= simp_penalty*Ta_DE[iVar];
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

          if (LockStrategy) omp_unset_lock(&UpdateLocks[indexNode[iNode]]);
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_MassMatrix(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Never record this method as the mass matrix is passive (but the mass residual is not). ---*/
  const bool ActiveTape = AD::TapeActive();
  AD::StopRecording();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Clear matrix before calculation. ---*/
    MassMatrix.SetValZero();

    for(auto color : ElemColoring) {

      /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
      SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];

        unsigned short iNode, jNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE_3D];

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

          if (LockStrategy) omp_set_lock(&UpdateLocks[indexNode[iNode]]);

          for (jNode = 0; jNode < nNodes; jNode++) {

            auto Mij = MassMatrix.GetBlock(indexNode[iNode], indexNode[jNode]);
            su2double Mab = simp_penalty * element->Get_Mab(iNode, jNode);

            for (iVar = 0; iVar < nVar; iVar++)
              Mij[iVar*(nVar+1)] += SU2_TYPE::GetValue(Mab);
          }

          if (LockStrategy) omp_unset_lock(&UpdateLocks[indexNode[iNode]]);
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

  if (ActiveTape) AD::StartRecording();

}

void CFEASolver::Compute_MassRes(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Clear vector before calculation. ---*/
  TimeRes.SetValZero();
  SU2_OMP_BARRIER

  for(auto color : ElemColoring) {

    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for(auto k = 0ul; k < color.size; ++k) {

      auto iElem = color.indices[k];

      unsigned short iNode, jNode, iDim, iVar;

      int thread = omp_get_thread_num();

      /*--- Convert VTK type to index in the element container. ---*/
      int EL_KIND;
      unsigned short nNodes;
      GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

      /*--- Each thread needs a dedicated element. ---*/
      CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

      /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
      unsigned long indexNode[MAXNNODE_3D];

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

      /*--- Add contributions of this element to the mass residual.
       *    Equiv. to a matrix vector product with the mass matrix. ---*/
      for (iNode = 0; iNode < nNodes; iNode++) {

        if (LockStrategy) omp_set_lock(&UpdateLocks[indexNode[iNode]]);

        for (jNode = 0; jNode < nNodes; jNode++) {

          su2double Mab = simp_penalty * element->Get_Mab(iNode, jNode);

          for (iVar = 0; iVar < nVar; iVar++)
            TimeRes(indexNode[iNode], iVar) += Mab * TimeRes_Aux(indexNode[jNode], iVar);
        }

        if (LockStrategy) omp_unset_lock(&UpdateLocks[indexNode[iNode]]);
      }

    } // end iElem loop

  } // end color loop

}

void CFEASolver::Compute_NodalStressRes(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  const bool prestretch_fem = config->GetPrestretch();

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();
  const su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Clear vector before calculation. ---*/
    LinSysRes.SetValZero();
    SU2_OMP_BARRIER

    for(auto color : ElemColoring) {

      /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
      SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];

        unsigned short iNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
        unsigned long indexNode[MAXNNODE_3D];

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
        int NUM_TERM = thread*MAX_TERMS + element_properties[iElem]->GetMat_Mod();

        numerics[NUM_TERM]->Compute_NodalStress_Term(element, config);

        for (iNode = 0; iNode < nNodes; iNode++) {
          if (LockStrategy) omp_set_lock(&UpdateLocks[indexNode[iNode]]);

          auto Ta = element->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes(indexNode[iNode], iVar) -= simp_penalty*Ta[iVar];

          if (LockStrategy) omp_unset_lock(&UpdateLocks[indexNode[iNode]]);
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_NodalStress(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  /*--- Never record this method as atm it is not differentiable. ---*/
  const bool ActiveTape = AD::TapeActive();
  AD::StopRecording();

  const bool prestretch_fem = config->GetPrestretch();

  const bool topology_mode = config->GetTopology_Optimization();
  const su2double simp_exponent = config->GetSIMP_Exponent();

  const unsigned short nStress = (nDim == 2) ? 3 : 6;

  su2double MaxVonMises_Stress = 0.0;

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Clear reactions. ---*/
    LinSysReact.SetValZero();

    /*--- Restart stress to avoid adding over results from previous time steps. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (unsigned short iStress = 0; iStress < nStress; iStress++) {
        nodes->SetStress_FEM(iPoint,iStress, 0.0);
      }
    }

    for(auto color : ElemColoring) {

      /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
      SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];

        unsigned short iNode, iDim, iVar, iStress;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
        unsigned long indexNode[MAXNNODE_3D];

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
        int NUM_TERM = thread*MAX_TERMS + element_properties[iElem]->GetMat_Mod();

        numerics[NUM_TERM]->Compute_Averaged_NodalStress(element, config);

        for (iNode = 0; iNode < nNodes; iNode++) {

          auto iPoint = indexNode[iNode];

          if (LockStrategy) omp_set_lock(&UpdateLocks[iPoint]);

          auto Ta = element->Get_Kt_a(iNode);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysReact(iPoint,iVar) += simp_penalty*Ta[iVar];

          /*--- Divide the nodal stress by the number of elements that will contribute to this point. ---*/
          su2double weight = simp_penalty / geometry->node[iPoint]->GetnElem();

          for (iStress = 0; iStress < nStress; iStress++)
            nodes->AddStress_FEM(iPoint,iStress, weight*element->Get_NodalStress(iNode,iStress));

          if (LockStrategy) omp_unset_lock(&UpdateLocks[iPoint]);
        }

      } // end iElem loop

    } // end color loop


    /*--- Compute the von Misses stress at each point, and the maximum for the domain. ---*/
    su2double maxVonMises = 0.0;

    SU2_OMP(for schedule(static,omp_chunk_size) nowait)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {

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

      maxVonMises = max(maxVonMises, VonMises_Stress);
    }
    SU2_OMP_CRITICAL
    MaxVonMises_Stress = max(MaxVonMises_Stress, maxVonMises);

  } // end SU2_OMP_PARALLEL

  su2double tmp = MaxVonMises_Stress;
  SU2_MPI::Allreduce(&tmp, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

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
                val_Reaction = LinSysReact(iPoint, iVar);
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
              TimeRes_Aux(iPoint,iVar) =
                a_dt[0]*nodes->GetSolution_time_n(iPoint,iVar) -      // a0*U(t)
                a_dt[0]*nodes->GetSolution(iPoint,iVar) +             // a0*U(t+dt)(k-1)
                a_dt[2]*nodes->GetSolution_Vel_time_n(iPoint,iVar) +  // a2*U'(t)
                a_dt[3]*nodes->GetSolution_Accel_time_n(iPoint,iVar); // a3*U''(t)
            }
          }

          /*--- Once computed, compute M*TimeRes_Aux ---*/
          Compute_MassRes(geometry, numerics, config);

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
                    val_Reaction = LinSysReact(iPoint, iVar) + TimeRes(iPoint, iVar);
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

  if (ActiveTape) AD::StartRecording();

}

void CFEASolver::Compute_DeadLoad(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- Clear integrated body forces before calculation. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      nodes->Clear_BodyForces_Res(iPoint);

    for(auto color : ElemColoring) {

      /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
      SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
      for(auto k = 0ul; k < color.size; ++k) {

        auto iElem = color.indices[k];

        unsigned short iNode, iDim, iVar;

        int thread = omp_get_thread_num();

        /*--- Convert VTK type to index in the element container. ---*/
        int EL_KIND;
        unsigned short nNodes;
        GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

        /*--- Each thread needs a dedicated element. ---*/
        CElement* element = element_container[FEA_TERM][EL_KIND+thread*MAX_FE_KINDS];

        /*--- For the number of nodes, get the coordinates and cache the point indices. ---*/
        unsigned long indexNode[MAXNNODE_3D];

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

          if (LockStrategy) omp_set_lock(&UpdateLocks[indexNode[iNode]]);

          auto Dead_Load = element->Get_FDL_a(iNode);

          su2double Aux_Dead_Load[MAXNVAR];
          for (iVar = 0; iVar < nVar; iVar++)
            Aux_Dead_Load[iVar] = simp_penalty*Dead_Load[iVar];

          nodes->Add_BodyForces_Res(indexNode[iNode], Aux_Dead_Load);

          if (LockStrategy) omp_unset_lock(&UpdateLocks[indexNode[iNode]]);
        }

      } // end iElem loop

    } // end color loop

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::Compute_IntegrationConstants(const CConfig *config) {

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


void CFEASolver::BC_Clamped(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  const bool dynamic = config->GetTime_Domain();
  const su2double zeros[MAXNVAR] = {0.0};

  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Set and enforce solution at current and previous time-step ---*/
    nodes->SetSolution(iPoint, zeros);

    if (dynamic) {
      nodes->SetSolution_Vel(iPoint, zeros);
      nodes->SetSolution_Accel(iPoint, zeros);
      nodes->Set_Solution_time_n(iPoint, zeros);
      nodes->SetSolution_Vel_time_n(iPoint, zeros);
      nodes->SetSolution_Accel_time_n(iPoint, zeros);
    }

    /*--- Set and enforce 0 solution for mesh deformation ---*/
    nodes->SetBound_Disp(iPoint, zeros);

    LinSysSol.SetBlock(iPoint, zeros);
    LinSysReact.SetBlock(iPoint, zeros);
    Jacobian.EnforceSolutionAtNode(iPoint, zeros, LinSysRes);

  }

}

void CFEASolver::BC_Clamped_Post(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  bool dynamic = config->GetTime_Domain();

  su2double zeros[MAXNVAR] = {0.0};

  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    nodes->SetSolution(iPoint, zeros);

    if (dynamic) {
      nodes->SetSolution_Vel(iPoint, zeros);
      nodes->SetSolution_Accel(iPoint, zeros);
    }

  }

}

void CFEASolver::BC_Sym_Plane(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  if (geometry->GetnElem_Bound(val_marker) == 0) return;

  const bool dynamic = config->GetTime_Domain();

  /*--- Determine axis of symmetry based on the normal of the first element in the marker. ---*/

  const su2double* nodeCoord[MAXNNODE_2D] = {nullptr};

  const bool quad = (geometry->bound[val_marker][0]->GetVTK_Type() == QUADRILATERAL);
  const unsigned short nNodes = quad? 4 : nDim;

  for (auto iNode = 0u; iNode < nNodes; iNode++) {
    auto iPoint = geometry->bound[val_marker][0]->GetNode(iNode);
    nodeCoord[iNode] = geometry->node[iPoint]->GetCoord();
  }

  su2double normal[MAXNDIM] = {0.0};

  switch (nNodes) {
    case 2: LineNormal(nodeCoord, normal); break;
    case 3: TriangleNormal(nodeCoord, normal); break;
    case 4: QuadrilateralNormal(nodeCoord, normal); break;
  }

  auto axis = 0u;
  for (auto iDim = 1u; iDim < MAXNDIM; ++iDim)
    axis = (fabs(normal[iDim]) > fabs(normal[axis]))? iDim : axis;

  if (fabs(normal[axis]) < 0.99*Norm(int(MAXNDIM),normal)) {
    SU2_MPI::Error("The structural solver only supports axis-aligned symmetry planes.",CURRENT_FUNCTION);
  }

  /*--- Impose zero displacement perpendicular to the symmetry plane. ---*/

  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Set and enforce solution at current and previous time-step ---*/
    nodes->SetSolution(iPoint, axis, 0.0);

    if (dynamic) {
      nodes->SetSolution_Vel(iPoint, axis, 0.0);
      nodes->SetSolution_Accel(iPoint, axis, 0.0);
      nodes->Set_Solution_time_n(iPoint, axis, 0.0);
      nodes->SetSolution_Vel_time_n(iPoint, axis, 0.0);
      nodes->SetSolution_Accel_time_n(iPoint, axis, 0.0);
    }

    /*--- Set and enforce 0 solution for mesh deformation ---*/
    nodes->SetBound_Disp(iPoint, axis, 0.0);

    LinSysSol(iPoint, axis) = 0.0;
    LinSysReact(iPoint, axis) = 0.0;
    Jacobian.EnforceSolutionAtDOF(iPoint, axis, su2double(0.0), LinSysRes);

  }

}

void CFEASolver::BC_DispDir(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  unsigned short iDim;

  auto TagBound = config->GetMarker_All_TagBound(val_marker);
  su2double DispDirVal = config->GetDisp_Dir_Value(TagBound);
  su2double DispDirMult = config->GetDisp_Dir_Multiplier(TagBound);
  const su2double *DispDirLocal = config->GetDisp_Dir(TagBound);
  su2double DispDirMod = Norm(nDim, DispDirLocal);

  su2double CurrentTime = config->GetCurrent_DynTime();
  su2double RampTime = config->GetRamp_Time();
  su2double ModAmpl = Compute_LoadCoefficient(CurrentTime, RampTime, config);

  su2double TotalDisp = ModAmpl * DispDirVal * DispDirMult / DispDirMod;

  su2double DispDir[MAXNVAR] = {0.0};
  for (iDim = 0; iDim < nDim; iDim++)
    DispDir[iDim] = TotalDisp * DispDirLocal[iDim];

  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index. ---*/
    auto iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- The solution is incremental so we need to
     *    subtract the current displacement. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      LinSysSol(iNode,iDim) = DispDir[iDim] - nodes->GetSolution(iNode,iDim);

    /*--- Enforce the solution. ---*/
    Jacobian.EnforceSolutionAtNode(iNode, LinSysSol.GetBlock(iNode), LinSysRes);
  }

}

void CFEASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                CConfig *config, CNumerics **numerics, unsigned short iMesh) {

  const bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);

  /*--- Compute stresses for monitoring and output. ---*/

  Compute_NodalStress(geometry, numerics, config);

  /*--- Compute the objective function to be able to monitor it. ---*/

  const auto kindObjFunc = config->GetKind_ObjFunc();

  if (((kindObjFunc == REFERENCE_GEOMETRY) || (kindObjFunc == REFERENCE_NODE)) &&
      ((config->GetDV_FEA() == YOUNG_MODULUS) || (config->GetDV_FEA() == DENSITY_VAL))) {

    Stiffness_Penalty(geometry, solver_container, numerics, config);
  }

  switch (kindObjFunc) {
    case REFERENCE_GEOMETRY: Compute_OFRefGeom(geometry, config); break;
    case REFERENCE_NODE:     Compute_OFRefNode(geometry, config); break;
    case VOLUME_FRACTION:    Compute_OFVolFrac(geometry, config); break;
    case TOPOL_DISCRETENESS: Compute_OFVolFrac(geometry, config); break;
    case TOPOL_COMPLIANCE:   Compute_OFCompliance(geometry, config); break;
  }

  if (nonlinear_analysis) {

    /*--- For nonlinear analysis we have 3 convergence criteria: ---*/
    /*--- UTOL = norm(Delta_U(k)): ABSOLUTE, norm of the incremental displacements ---*/
    /*--- RTOL = norm(Residual(k): ABSOLUTE, norm of the residual (T-F) ---*/
    /*--- ETOL = Delta_U(k) * Residual(k): ABSOLUTE, energy norm ---*/

    SU2_OMP_PARALLEL
    {
    su2double utol = LinSysSol.norm();
    su2double rtol = LinSysRes.norm();
    su2double etol = LinSysSol.dot(LinSysRes);

    SU2_OMP_MASTER
    {
      Conv_Check[0] = utol;
      Conv_Check[1] = rtol;
      Conv_Check[2] = etol;
    }
    } // end parallel
  }
  else {

    /*--- If the problem is linear, the only check we do is the RMS of the residuals. ---*/
    /*---  Compute the residual Ax-f ---*/

#ifndef CODI_FORWARD_TYPE
    CSysVector<passivedouble> LinSysAux(nPoint, nPointDomain, nVar, nullptr);
#else
    CSysVector<su2double> LinSysAux(nPoint, nPointDomain, nVar, nullptr);
#endif

    SU2_OMP_PARALLEL
    {
#ifndef CODI_REVERSE_TYPE
    Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
#else
    /*---  We need temporaries to interface with the passive matrix. ---*/
    CSysVector<passivedouble> sol, res;
    sol.PassiveCopy(LinSysSol);
    res.PassiveCopy(LinSysRes);
    Jacobian.ComputeResidual(sol, res, LinSysAux);
#endif

    /*--- Set maximum residual to zero. ---*/

    SU2_OMP_MASTER
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
      SetRes_RMS(iVar, 0.0);
      SetRes_Max(iVar, 0.0, 0);
    }

    /*--- Compute the residual. ---*/

    su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
    const su2double* coordMax[MAXNVAR] = {nullptr};
    unsigned long idxMax[MAXNVAR] = {0};

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      for (auto iVar = 0ul; iVar < nVar; iVar++) {
        su2double Res = fabs(LinSysAux(iPoint, iVar));
        resRMS[iVar] += Res*Res;
        if (Res > resMax[iVar]) {
          resMax[iVar] = Res;
          idxMax[iVar] = iPoint;
          coordMax[iVar] = geometry->node[iPoint]->GetCoord();
        }
      }
    }
    SU2_OMP_CRITICAL
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
      AddRes_RMS(iVar, resRMS[iVar]);
      AddRes_Max(iVar, resMax[iVar], geometry->node[idxMax[iVar]]->GetGlobalIndex(), coordMax[iVar]);
    }
    SU2_OMP_BARRIER

    /*--- Compute the root mean square residual. ---*/
    SU2_OMP_MASTER
    SetResidual_RMS(geometry, config);

    } // end SU2_OMP_PARALLEL

  }

}

void CFEASolver::BC_Normal_Load(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  /*--- Determine whether the load conditions are applied in the reference or in the current configuration. ---*/

  const bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);

  /*--- Retrieve the normal pressure and the application conditions for the considered boundary. ---*/

  su2double CurrentTime = config->GetCurrent_DynTime();
  su2double Ramp_Time = config->GetRamp_Time();
  su2double ModAmpl = Compute_LoadCoefficient(CurrentTime, Ramp_Time, config);

  su2double NormalLoad = config->GetLoad_Value(config->GetMarker_All_TagBound(val_marker));
  const su2double TotalLoad = ModAmpl * NormalLoad;

  /*--- Continue only if there is a load applied, to reduce computational cost. ---*/

  if (TotalLoad == 0.0) return;


  for (unsigned long iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    unsigned short iNode, iDim;
    unsigned long indexNode[MAXNNODE_2D] = {0};
    su2double nodeCoord_ref[MAXNNODE_2D][MAXNDIM] = {{0.0}};
    su2double nodeCoord_curr[MAXNNODE_2D][MAXNDIM] = {{0.0}};

    /*--- Identify the kind of boundary element. ---*/

    bool quad = (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL);
    unsigned short nNodes = quad? 4 : nDim;

    /*--- Retrieve the boundary reference and current coordinates. ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        nodeCoord_ref[iNode][iDim] = val_Coord;

        if (nonlinear_analysis) {
          val_Coord += nodes->GetSolution(indexNode[iNode],iDim);
          nodeCoord_curr[iNode][iDim] = val_Coord;
        }
      }
    }

    /*--- Compute area vectors in reference and current configurations. ---*/

    su2double normal_ref[MAXNDIM] = {0.0};
    su2double normal_curr[MAXNDIM] = {0.0};

    switch (nNodes) {
      case 2: LineNormal(nodeCoord_ref, normal_ref); break;
      case 3: TriangleNormal(nodeCoord_ref, normal_ref); break;
      case 4: QuadrilateralNormal(nodeCoord_ref, normal_ref); break;
    }

    if (nonlinear_analysis) {
      switch (nNodes) {
        case 2: LineNormal(nodeCoord_curr, normal_curr); break;
        case 3: TriangleNormal(nodeCoord_curr, normal_curr); break;
        case 4: QuadrilateralNormal(nodeCoord_curr, normal_curr); break;
      }
    }

    /*--- Use a reference normal from one of the points to decide if computed normal needs to be flipped. ---*/

    auto reference_vertex = geometry->node[indexNode[0]]->GetVertex(val_marker);
    const su2double* reference_normal = geometry->vertex[val_marker][reference_vertex]->GetNormal();

    su2double dot = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dot += normal_ref[iDim] * reference_normal[iDim];

    su2double sign = (dot < 0.0)? -1.0 : 1.0;

    /*--- Compute the load vector (overwrites one of the normals). ---*/

    su2double* load = nonlinear_analysis? normal_curr : normal_ref;

    for (iDim = 0; iDim < nDim; iDim++)
      load[iDim] *= sign * TotalLoad / su2double(nNodes);

    /*--- Update surface load for each node of the boundary element. ---*/

    for (iNode = 0; iNode < nNodes; iNode++)
      nodes->Add_SurfaceLoad_Res(indexNode[iNode], load);
  }

}

void CFEASolver::BC_Dir_Load(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  auto TagBound = config->GetMarker_All_TagBound(val_marker);
  su2double LoadDirVal = config->GetLoad_Dir_Value(TagBound);
  su2double LoadDirMult = config->GetLoad_Dir_Multiplier(TagBound);
  const su2double* Load_Dir_Local = config->GetLoad_Dir(TagBound);

  /*--- Compute the norm of the vector that was passed in the config file. ---*/
  su2double LoadNorm = Norm(nDim, Load_Dir_Local);

  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double Ramp_Time = config->GetRamp_Time();
  su2double ModAmpl = Compute_LoadCoefficient(CurrentTime, Ramp_Time, config);

  const su2double TotalLoad = ModAmpl * LoadDirVal * LoadDirMult / LoadNorm;


  for (unsigned long iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    unsigned short iNode, iDim;
    unsigned long indexNode[MAXNNODE_2D] = {0};

    const su2double* nodeCoord[MAXNNODE_2D] = {nullptr};

    /*--- Identify the kind of boundary element. ---*/

    bool quad = (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL);
    unsigned short nNodes = quad? 4 : nDim;

    /*--- Retrieve the boundary reference coordinates. ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);
      nodeCoord[iNode] = geometry->node[indexNode[iNode]]->GetCoord();
    }

    /*--- Compute area of the boundary element. ---*/

    su2double normal[MAXNDIM] = {0.0};

    switch (nNodes) {
      case 2: LineNormal(nodeCoord, normal); break;
      case 3: TriangleNormal(nodeCoord, normal); break;
      case 4: QuadrilateralNormal(nodeCoord, normal); break;
    }

    su2double area = Norm(int(MAXNDIM),normal);

    /*--- Compute load vector and update surface load for each node of the boundary element. ---*/

    su2double* load = normal;

    for (iDim = 0; iDim < nDim; ++iDim)
      load[iDim] = Load_Dir_Local[iDim] * TotalLoad * area / su2double(nNodes);

    for (iNode = 0; iNode < nNodes; iNode++)
      nodes->Add_SurfaceLoad_Res(indexNode[iNode], load);
  }

}

void CFEASolver::BC_Damper(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker) {

  const su2double dampConst = config->GetDamper_Constant(config->GetMarker_All_TagBound(val_marker));

  for (auto iElem = 0ul; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    unsigned short iNode, iDim;
    unsigned long indexNode[MAXNNODE_2D] = {0};

    su2double nodeCoord[MAXNNODE_2D][MAXNDIM] = {{0.0}};

    bool quad = (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL);
    unsigned short nNodes = quad? 4 : nDim;

    /*--- Retrieve the boundary current coordinates. ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {

      auto iPoint = geometry->bound[val_marker][iElem]->GetNode(iNode);
      indexNode[iNode] = iPoint;

      for (iDim = 0; iDim < nVar; iDim++)
        nodeCoord[iNode][iDim] = geometry->node[iPoint]->GetCoord(iDim) + nodes->GetSolution(iPoint,iDim);
    }

    /*--- Compute the area of the surface element. ---*/

    su2double normal[MAXNDIM] = {0.0};

    switch (nNodes) {
      case 2: LineNormal(nodeCoord, normal); break;
      case 3: TriangleNormal(nodeCoord, normal); break;
      case 4: QuadrilateralNormal(nodeCoord, normal); break;
    }

    su2double area = Norm(int(MAXNDIM),normal);

    /*--- Compute damping forces. ---*/

    su2double dampCoeff = -1.0 * area * dampConst / su2double(nNodes);

    for(iNode = 0; iNode < nNodes; ++iNode) {

      auto iPoint = indexNode[iNode];

      /*--- Writing over the normal. --*/
      su2double* force = normal;
      for (iDim = 0; iDim < nVar; iDim++)
        force[iDim] = dampCoeff * nodes->GetSolution_Vel(iPoint, iDim);

      nodes->Add_SurfaceLoad_Res(iPoint, force);
    }

  }

}

void CFEASolver::BC_Deforming(CGeometry *geometry, CNumerics *numerics, const CConfig *config, unsigned short val_marker){

  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    auto iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Retrieve the boundary displacement ---*/
    su2double Disp[MAXNVAR] = {0.0};
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Disp[iDim] = nodes->GetBound_Disp(iNode,iDim);

    /*--- Set and enforce solution ---*/
    LinSysSol.SetBlock(iNode, Disp);
    Jacobian.EnforceSolutionAtNode(iNode, Disp, LinSysRes);

  }

}

void CFEASolver::Integrate_FSI_Loads(CGeometry *geometry, const CConfig *config) {

  /*--- The conservative approach transfers forces directly, no integration needed. ---*/
  if (config->GetConservativeInterpolation()) return;

  unordered_map<unsigned long, su2double> vertexArea;
  unordered_set<short> processedMarkers({-1});

  /*--- Compute current area associated with each vertex. ---*/

  for (auto iMarkerInt = 0; iMarkerInt < config->GetMarker_n_ZoneInterface()/2; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair. ---*/
    const auto iMarker = config->FindInterfaceMarker(iMarkerInt);
    /*--- The current mpi rank may not have this marker, or it may have been processed already. ---*/
    if (processedMarkers.count(iMarker) > 0) continue;
    processedMarkers.insert(iMarker);

    for (auto iElem = 0u; iElem < geometry->GetnElem_Bound(iMarker); ++iElem) {
      /*--- Define the boundary element. ---*/
      unsigned long nodeList[MAXNNODE_2D] = {0};
      su2double coords[MAXNNODE_2D][MAXNDIM] = {{0.0}};
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      auto nNode = quad? 4u : nDim;

      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        nodeList[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          coords[iNode][iDim] = geometry->node[nodeList[iNode]]->GetCoord(iDim)+
                                nodes->GetSolution(nodeList[iNode],iDim);
      }

      /*--- Compute the area contribution to each node. ---*/
      su2double normal[MAXNDIM] = {0.0};

      switch (nNode) {
        case 2: LineNormal(coords, normal); break;
        case 3: TriangleNormal(coords, normal); break;
        case 4: QuadrilateralNormal(coords, normal); break;
      }

      su2double area = Norm(int(MAXNDIM),normal) / nNode;

      /*--- Update area of nodes. ---*/
      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        auto iPoint = nodeList[iNode];
        if (vertexArea.count(iPoint) == 0)
          vertexArea[iPoint] = area;
        else
          vertexArea[iPoint] += area;
      }
    }
  }

  /*--- Integrate tractions. ---*/

  for (auto it = vertexArea.begin(); it != vertexArea.end(); ++it) {
    auto iPoint = it->first;
    su2double area = it->second;
    su2double force[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; ++iDim)
      force[iDim] = nodes->Get_FlowTraction(iPoint,iDim)*area;
    nodes->Set_FlowTraction(iPoint, force);
  }

}

su2double CFEASolver::Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, const CConfig *config){

  su2double LoadCoeff = 1.0;

  bool Ramp_Load = config->GetRamp_Load();
  bool Sine_Load = config->GetSine_Load();
  bool Ramp_And_Release = config->GetRampAndRelease_Load();

  bool restart = config->GetRestart(); // Restart analysis
  bool fsi = config->GetFSI_Simulation();
  bool stat_fsi = !config->GetTime_Domain();

  /*--- This offset introduces the ramp load in dynamic cases starting from the restart point. ---*/
  bool offset = (restart && fsi && (!stat_fsi));
  su2double DeltaT = config->GetDelta_DynTime();
  su2double OffsetTime = offset? DeltaT * (config->GetRestart_Iter()-1) : su2double(0.0);

  /*--- Polynomial functions from https://en.wikipedia.org/wiki/Smoothstep ---*/

  if ((Ramp_Load) && (RampTime > 0.0)) {

    su2double TransferTime = (CurrentTime - OffsetTime) / RampTime;

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
    su2double SineAmp   = config->GetLoad_Sine()[0];
    su2double SineFreq  = config->GetLoad_Sine()[1];
    su2double SinePhase = config->GetLoad_Sine()[2];

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

void CFEASolver::ImplicitNewmark_Iteration(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

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
        LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_SurfaceLoad_Res(iPoint,iVar);
      }

      /*--- Body forces contribution (dead load). ---*/

      if (body_forces) {
        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_BodyForces_Res(iPoint,iVar);
        }
      }

      /*--- FSI contribution (flow loads). ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_FlowTraction(iPoint,iVar);
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
        Jacobian.MatrixMatrixAddition(SU2_TYPE::GetValue(a_dt[0]), MassMatrix);
      }

      /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          TimeRes_Aux(iPoint,iVar) =
            a_dt[0]*nodes->GetSolution_time_n(iPoint,iVar) -      // a0*U(t)
            a_dt[0]*nodes->GetSolution(iPoint,iVar) +             // a0*U(t+dt)(k-1)
            a_dt[2]*nodes->GetSolution_Vel_time_n(iPoint,iVar) +  // a2*U'(t)
            a_dt[3]*nodes->GetSolution_Accel_time_n(iPoint,iVar); // a3*U''(t)
        }
      }

      /*--- Add M*TimeRes_Aux to the residual. ---*/
      Compute_MassRes(geometry, numerics, config);
      LinSysRes += TimeRes;
    }

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::ImplicitNewmark_Update(CGeometry *geometry, CConfig *config) {

  const bool dynamic = (config->GetTime_Domain());

  SU2_OMP_PARALLEL
  {
    unsigned long iPoint;
    unsigned short iVar;

    /*--- Update solution. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Displacement component of the solution. ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        nodes->Add_DeltaSolution(iPoint, iVar, LinSysSol(iPoint,iVar));
    }

    if (dynamic) {
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {

          /*--- Acceleration component of the solution. ---*/
          /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

          su2double sol = a_dt[0]*(nodes->GetSolution(iPoint,iVar) -
                                   nodes->GetSolution_time_n(iPoint,iVar)) -
                          a_dt[2]* nodes->GetSolution_Vel_time_n(iPoint,iVar) -
                          a_dt[3]* nodes->GetSolution_Accel_time_n(iPoint,iVar);

          nodes->SetSolution_Accel(iPoint, iVar, sol);

          /*--- Velocity component of the solution. ---*/
          /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

          sol = nodes->GetSolution_Vel_time_n(iPoint,iVar)+
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

void CFEASolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CConfig *config) {

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

          /*--- Velocity component of the solution ---*/
          /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

          sol = nodes->GetSolution_Vel_time_n(iPoint,iVar)+
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


void CFEASolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CNumerics **numerics, const CConfig *config) {

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
          LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_SurfaceLoad_Res(iPoint,iVar);
        }

        /*--- Body forces contribution (dead load). ---*/

        if (body_forces) {
          for (iVar = 0; iVar < nVar; iVar++) {
            LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_BodyForces_Res(iPoint,iVar);
          }
        }

        /*--- FSI contribution (flow loads). ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_FlowTraction(iPoint,iVar);
        }

      }
    }

    /*--- Loads for dynamic problems. ---*/

    if (dynamic) {

      /*--- Add the mass matrix contribution to the Jacobian. ---*/

      /*--- See notes on logic in ImplicitNewmark_Iteration(). ---*/
      if ((nonlinear_analysis && (newton_raphson || first_iter)) || linear_analysis) {
        Jacobian.MatrixMatrixAddition(SU2_TYPE::GetValue(a_dt[0]), MassMatrix);
      }

      /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          TimeRes_Aux(iPoint,iVar) =
            a_dt[0]*nodes->GetSolution_time_n(iPoint,iVar) -      // a0*U(t)
            a_dt[0]*nodes->GetSolution(iPoint,iVar) +             // a0*U(t+dt)(k-1)
            a_dt[2]*nodes->GetSolution_Vel_time_n(iPoint,iVar) +  // a2*U'(t)
            a_dt[3]*nodes->GetSolution_Accel_time_n(iPoint,iVar); // a3*U''(t)
        }
      }

      /*--- Add M*TimeRes_Aux to the residual. ---*/
      Compute_MassRes(geometry, numerics, config);
      LinSysRes += TimeRes;
      SU2_OMP_BARRIER

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {

        /*--- External surface load contribution ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes(iPoint,iVar) += loadIncr * ( (1-alpha_f) * nodes->Get_SurfaceLoad_Res(iPoint,iVar) +
                                                    alpha_f  * nodes->Get_SurfaceLoad_Res_n(iPoint,iVar) );
        }

        /*--- Add the contribution to the residual due to body forces.
         *--- It is constant over time, so it's not necessary to distribute it. ---*/

        if (body_forces) {
          for (iVar = 0; iVar < nVar; iVar++) {
            LinSysRes(iPoint,iVar) += loadIncr * nodes->Get_BodyForces_Res(iPoint,iVar);
          }
        }

        /*--- Add FSI contribution. ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes(iPoint,iVar) += loadIncr * ( (1-alpha_f) * nodes->Get_FlowTraction(iPoint,iVar) +
                                                    alpha_f  * nodes->Get_FlowTraction_n(iPoint,iVar) );
        }
      }
    }

  } // end SU2_OMP_PARALLEL

}

void CFEASolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CConfig *config) {

  /*--- Update displacement components of the solution. ---*/

  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size))
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      nodes->Add_DeltaSolution(iPoint, iVar, LinSysSol(iPoint,iVar));

  /*--- Perform the MPI communication of the solution, displacements only. ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

}

void CFEASolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CConfig *config) {

  const su2double alpha_f = config->Get_Int_Coeffs(2);
  const su2double alpha_m = config->Get_Int_Coeffs(3);

  /*--- Compute solution at t_n+1, and update velocities and accelerations ---*/

  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size))
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

      /*--- Velocity component of the solution ---*/
      /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

      sol = nodes->GetSolution_Vel_time_n(iPoint,iVar)+
            a_dt[6]* nodes->GetSolution_Accel_time_n(iPoint,iVar) +
            a_dt[7]* nodes->GetSolution_Accel(iPoint,iVar);

      nodes->SetSolution_Vel(iPoint, iVar, sol);
    }

  }

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION_FEA);
  CompleteComms(geometry, config, SOLUTION_FEA);

}

void CFEASolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, const CConfig *config) {

  /*--- Set the load conditions of the time step n+1 as the load conditions for time step n ---*/
  nodes->Set_SurfaceLoad_Res_n();
  nodes->Set_FlowTraction_n();

}

void CFEASolver::Solve_System(CGeometry *geometry, CConfig *config) {

  /*--- Enforce solution at some halo points possibly not covered by essential BC markers. ---*/

  Jacobian.InitiateComms(LinSysSol, geometry, config, SOLUTION_MATRIX);
  Jacobian.CompleteComms(LinSysSol, geometry, config, SOLUTION_MATRIX);

  for (auto iPoint : ExtraVerticesToEliminate) {
    Jacobian.EnforceSolutionAtNode(iPoint, LinSysSol.GetBlock(iPoint), LinSysRes);
  }

  SU2_OMP_PARALLEL
  {
  /*--- This is required for the discrete adjoint. ---*/
  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto i = nPointDomain*nVar; i < nPoint*nVar; ++i) LinSysRes[i] = 0.0;

  /*--- Solve or smooth the linear system. ---*/

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  SU2_OMP_MASTER
  {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  //SU2_OMP_BARRIER
  } // end SU2_OMP_PARALLEL

}


void CFEASolver::PredictStruct_Displacement(CGeometry *geometry, CConfig *config) {

  const unsigned short predOrder = config->GetPredictorOrder();
  const su2double Delta_t = config->GetDelta_DynTime();

  if(predOrder > 2 && rank == MASTER_NODE)
    cout << "Higher order predictor not implemented. Solving with order 0." << endl;

  /*--- To nPointDomain: we need to communicate the predicted solution after setting it. ---*/
  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size))
  for (unsigned long iPoint=0; iPoint < nPointDomain; iPoint++) {

    unsigned short iDim;

    switch (predOrder) {
      case 1: {
        const su2double* solDisp = nodes->GetSolution(iPoint);
        const su2double* solVel = nodes->GetSolution_Vel(iPoint);
        su2double* valPred = nodes->GetSolution_Pred(iPoint);

        for (iDim=0; iDim < nDim; iDim++) {
          valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
        }
      } break;

      case 2: {
        const su2double* solDisp = nodes->GetSolution(iPoint);
        const su2double* solVel = nodes->GetSolution_Vel(iPoint);
        const su2double* solVel_tn = nodes->GetSolution_Vel_time_n(iPoint);
        su2double* valPred = nodes->GetSolution_Pred(iPoint);

        for (iDim=0; iDim < nDim; iDim++) {
          valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
        }
      } break;

      default: {
        nodes->SetSolution_Pred(iPoint);
      } break;
    }

  }

  InitiateComms(geometry, config, SOLUTION_PRED);
  CompleteComms(geometry, config, SOLUTION_PRED);

}

void CFEASolver::ComputeAitken_Coefficient(CGeometry *geometry, CConfig *config, unsigned long iOuterIter) {

  unsigned long iPoint, iDim;
  su2double rbuf_numAitk = 0, sbuf_numAitk = 0;
  su2double rbuf_denAitk = 0, sbuf_denAitk = 0;

  const su2double *dispPred = nullptr;
  const su2double *dispCalc = nullptr;
  const su2double *dispPred_Old = nullptr;
  const su2double *dispCalc_Old = nullptr;
  su2double deltaU[MAXNVAR] = {0.0}, deltaU_p1[MAXNVAR] = {0.0};
  su2double delta_deltaU[MAXNVAR] = {0.0};
  su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn_Min, WAitkDyn;

  unsigned short RelaxMethod_FSI = config->GetRelaxation_Method_FSI();

  /*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/

  if (RelaxMethod_FSI == NO_RELAXATION) {

    SetWAitken_Dyn(1.0);

  }
  else if (RelaxMethod_FSI == FIXED_PARAMETER) {

    SetWAitken_Dyn(config->GetAitkenStatRelax());

  }
  else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {

    if (iOuterIter == 0) {

      WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
      WAitkDyn_Max = config->GetAitkenDynMaxInit();
      WAitkDyn_Min = config->GetAitkenDynMinInit();

      WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);
      WAitkDyn = max(WAitkDyn, WAitkDyn_Min);

      SetWAitken_Dyn(WAitkDyn);

    }
    else {
      // To nPointDomain; we need to communicate the values
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        dispPred     = nodes->GetSolution_Pred(iPoint);
        dispPred_Old = nodes->GetSolution_Pred_Old(iPoint);
        dispCalc     = nodes->GetSolution(iPoint);
        dispCalc_Old = nodes->GetSolution_Old(iPoint);

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

void CFEASolver::SetAitken_Relaxation(CGeometry *geometry, CConfig *config) {

  const su2double WAitken = GetWAitken_Dyn();

  // To nPointDomain; we need to communicate the solutions (predicted, old and old predicted) after this routine
  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size))
  for (unsigned long iPoint=0; iPoint < nPointDomain; iPoint++) {

    /*--- Retrieve pointers to the predicted and calculated solutions ---*/
    su2double* dispPred = nodes->GetSolution_Pred(iPoint);
    const su2double* dispCalc = nodes->GetSolution(iPoint);

    /*--- Set predicted solution as the old predicted solution ---*/
    nodes->SetSolution_Pred_Old(iPoint);

    /*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
    nodes->SetSolution_Old(iPoint, dispCalc);

    /*--- Apply the Aitken relaxation ---*/
    for (unsigned short iDim=0; iDim < nDim; iDim++) {
      dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
    }
  }

  InitiateComms(geometry, config, SOLUTION_PRED_OLD);
  CompleteComms(geometry, config, SOLUTION_PRED_OLD);

}

void CFEASolver::OutputForwardModeGradient(const CConfig *config, bool newFile,
                                           su2double fun, su2double fun_avg,
                                           su2double der, su2double der_avg) const {
  if (rank != MASTER_NODE) return;

  bool dynamic = config->GetTime_Domain();

  string fileSuffix, varName;

  switch (config->GetDirectDiff()) {
    case D_YOUNG:
      fileSuffix = "E";
      varName = "Young's modulus";
      break;
    case D_POISSON:
      fileSuffix = "Nu";
      varName = "Poisson's ratio";
      break;
    case D_RHO:
      fileSuffix = "Rho";
      varName = "structural density";
      break;
    case D_RHO_DL:
      fileSuffix = "Rho_DL";
      varName = "dead weight";
      break;
    case D_EFIELD:
      fileSuffix = "EField";
      varName = "electric field";
      break;
    case D_MACH:
      fileSuffix = "Mach";
      varName = "Mach number";
      break;
    case D_PRESSURE:
      fileSuffix = "Pressure";
      varName = "freestream pressure";
      break;
    default:
      return;
      break;
  }

  ofstream myfile_res;

  if (newFile) {
    myfile_res.open(string("Output_Direct_Diff_")+fileSuffix+string(".txt"));

    myfile_res << "Objective Function" << "\t";
    if (dynamic)
      myfile_res << "O. Function Averaged" << "\t";
    myfile_res << "Sensitivity Local" << "\t";
    myfile_res << "Sensitivity Averaged" << endl;
    return;
  }

  myfile_res.open(string("Output_Direct_Diff_")+fileSuffix+string(".txt"), ios::app);
  myfile_res.precision(15);

  myfile_res << scientific << fun << "\t";
  if (dynamic)
    myfile_res << scientific << fun_avg << "\t";
  myfile_res << scientific << der << "\t";
  myfile_res << scientific << der_avg << endl;

  cout << "Objective function: " << fun << ". Global derivative of the "
       << varName << ": " << Total_ForwardGradient << "." << endl;

}

void CFEASolver::Compute_OFRefGeom(CGeometry *geometry, const CConfig *config){

  bool fsi = config->GetFSI_Simulation();
  unsigned long TimeIter = config->GetTimeIter();

  su2double objective_function = 0.0;

  SU2_OMP_PARALLEL
  {
  su2double obj_fun_local = 0.0;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {

      /*--- Retrieve the value of the reference geometry ---*/
      su2double reference_geometry = nodes->GetReference_Geometry(iPoint,iVar);

      /*--- Retrieve the value of the current solution ---*/
      su2double current_solution = nodes->GetSolution(iPoint,iVar);

      /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
      obj_fun_local += pow(current_solution - reference_geometry, 2);
    }
  }
  atomicAdd(obj_fun_local, objective_function);
  }

  SU2_MPI::Allreduce(&objective_function, &Total_OFRefGeom, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Total_OFRefGeom *= config->GetRefGeom_Penalty() / geometry->GetGlobal_nPointDomain();
  Total_OFRefGeom += PenaltyValue;

  Global_OFRefGeom += Total_OFRefGeom;

  /*--- To be accessible from the output. ---*/
  Total_OFCombo = Total_OFRefGeom;

  /// TODO: Temporary output files for the direct mode.

  if ((rank == MASTER_NODE) && (config->GetDirectDiff() != NO_DERIVATIVE)) {

    /*--- Forward mode AD results. ---*/

    su2double local_forward_gradient = SU2_TYPE::GetDerivative(Total_OFRefGeom);
    su2double objective_function_averaged = Global_OFRefGeom / (TimeIter + 1.0 + EPS);

    if (fsi) Total_ForwardGradient  = local_forward_gradient;
    else     Total_ForwardGradient += local_forward_gradient;

    su2double averaged_gradient = Total_ForwardGradient / (TimeIter + 1.0);

    OutputForwardModeGradient(config, false, Total_OFRefGeom, objective_function_averaged,
                              local_forward_gradient, averaged_gradient);
  }

}

void CFEASolver::Compute_OFRefNode(CGeometry *geometry, const CConfig *config){

  bool fsi = config->GetFSI_Simulation();
  unsigned long TimeIter = config->GetTimeIter();

  su2double dist[MAXNVAR] = {0.0}, dist_reduce[MAXNVAR];

  /*--- Convert global point index to local. ---*/
  long iPoint = geometry->GetGlobal_to_Local_Point(config->GetRefNode_ID());

  if (iPoint >= 0) {
    if (geometry->node[iPoint]->GetDomain()) {
      for (unsigned short iVar = 0; iVar < nVar; ++iVar)
        dist[iVar] = nodes->GetSolution(iPoint,iVar) - config->GetRefNode_Displacement(iVar);
    }
  }

  SU2_MPI::Allreduce(dist, dist_reduce, MAXNVAR, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Total_OFRefNode = config->GetRefNode_Penalty() * Norm(int(MAXNVAR),dist_reduce) + PenaltyValue;

  Global_OFRefNode += Total_OFRefNode;

  /*--- To be accessible from the output. ---*/
  Total_OFCombo = Total_OFRefNode;

  /// TODO: Temporary output files for the direct mode.

  if ((rank == MASTER_NODE) && (config->GetDirectDiff() != NO_DERIVATIVE)) {

    /*--- Forward mode AD results. ---*/

    su2double local_forward_gradient = SU2_TYPE::GetDerivative(Total_OFRefNode);
    su2double objective_function_averaged = Global_OFRefNode / (TimeIter + 1.0 + EPS);

    if (fsi) Total_ForwardGradient  = local_forward_gradient;
    else     Total_ForwardGradient += local_forward_gradient;

    su2double averaged_gradient = Total_ForwardGradient / (TimeIter + 1.0);

    OutputForwardModeGradient(config, false, Total_OFRefNode, objective_function_averaged,
                              local_forward_gradient, averaged_gradient);
  }

}

void CFEASolver::Compute_OFVolFrac(CGeometry *geometry, const CConfig *config)
{
  /*--- Perform a volume average of the physical density of the elements for topology optimization ---*/

  su2double total_volume = 0.0, integral = 0.0, discreteness = 0.0;

  SU2_OMP_PARALLEL
  {
  su2double tot_vol_loc = 0.0, integral_loc = 0.0, discrete_loc = 0.0;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iElem = 0; iElem < nElement; ++iElem) {
    /*--- count only elements that belong to the partition ---*/
    if (geometry->node[geometry->elem[iElem]->GetNode(0)]->GetDomain()) {
      su2double volume = geometry->elem[iElem]->GetVolume();
      su2double rho = element_properties[iElem]->GetPhysicalDensity();
      tot_vol_loc += volume;
      integral_loc += volume*rho;
      discrete_loc += volume*4.0*rho*(1.0-rho);
    }
  }
  atomicAdd(tot_vol_loc, total_volume);
  atomicAdd(integral_loc, integral);
  atomicAdd(discrete_loc, discreteness);
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

  /*--- To be accessible from the output. ---*/
  Total_OFCombo = Total_OFVolFrac;

}

void CFEASolver::Compute_OFCompliance(CGeometry *geometry, const CConfig *config)
{
  /*--- Types of loads to consider ---*/
  const bool body_forces = config->GetDeadLoad();

  /*--- If the loads are being applied incrementaly ---*/
  const bool incremental_load = config->GetIncrementalLoad();

  /*--- Computation (compliance = sum(f dot u) ) ---*/
  /*--- Cannot be computed as u^T K u as the problem may not be linear ---*/

  su2double compliance = 0.0;

  SU2_OMP_PARALLEL
  {
  su2double comp_local = 0.0;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    unsigned short iVar;
    su2double nodalForce[MAXNVAR];

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
      comp_local += nodalForce[iVar]*nodes->GetSolution(iPoint,iVar);
  }
  atomicAdd(comp_local, compliance);
  }

  SU2_MPI::Allreduce(&compliance, &Total_OFCompliance, 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  /*--- To be accessible from the output. ---*/
  Total_OFCombo = Total_OFCompliance;

}

void CFEASolver::Stiffness_Penalty(CGeometry *geometry, CSolver **solver, CNumerics **numerics, CConfig *config){

  if (config->GetTotalDV_Penalty() == 0.0) {
    /*--- No need to go into expensive computations. ---*/
    PenaltyValue = 0.0;
    return;
  }

  su2double weightedValue = 0.0;
  su2double weightedValue_reduce = 0.0;
  su2double totalVolume = 0.0;
  su2double totalVolume_reduce = 0.0;

  /*--- Loop over the elements in the domain. ---*/
  SU2_OMP_PARALLEL
  {
  su2double weighted_loc = 0.0, totalVol_loc = 0.0;

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iElem = 0; iElem < nElement; iElem++) {

    int thread = omp_get_thread_num();

    int EL_KIND;
    unsigned short iNode, nNodes, iDim;
    unsigned long indexNode[MAXNNODE_3D];

    GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

    CElement* element = element_container[FEA_TERM][EL_KIND + thread*MAX_FE_KINDS];

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        su2double val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element->SetRef_Coord(iNode, iDim, val_Coord);
      }
    }

    // Avoid double-counting elements:
    // Only add the value if the first node is in the domain
    if (geometry->node[indexNode[0]]->GetDomain()) {

      // Compute the area/volume of the element
      su2double elementVolume;

      if (nDim == 2)
        elementVolume = element->ComputeArea();
      else
        elementVolume = element->ComputeVolume();

      // Compute the total volume
      totalVol_loc += elementVolume;

      // Retrieve the value of the design variable
      su2double dvValue = numerics[FEA_TERM]->Get_DV_Val(element_properties[iElem]->GetDV());

      // Add the weighted sum of the value of the design variable
      weighted_loc += dvValue * elementVolume;

    }
  }
  atomicAdd(totalVol_loc, totalVolume);
  atomicAdd(weighted_loc, weightedValue);
  }

  // Reduce value across processors for parallelization

  SU2_MPI::Allreduce(&weightedValue, &weightedValue_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&totalVolume, &totalVolume_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  su2double ratio = 1.0 - weightedValue_reduce/totalVolume_reduce;

  PenaltyValue = config->GetTotalDV_Penalty() * ratio * ratio;

}

void CFEASolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  const bool dynamic = (config->GetTime_Domain());
  const bool fluid_structure = config->GetFSI_Simulation();
  const bool discrete_adjoint = config->GetDiscrete_Adjoint();

  /*--- Skip coordinates ---*/

  const auto skipVars = geometry[MESH_0]->GetnDim();

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  string filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  unsigned long iPoint_Global, counter = 0;

  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local >= 0) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      const auto index = counter*Restart_Vars[1] + skipVars;
      const passivedouble* Sol = &Restart_Data[index];

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        nodes->SetSolution(iPoint_Local, iVar, Sol[iVar]);
        if (dynamic) {
          nodes->Set_Solution_time_n(iPoint_Local, iVar, Sol[iVar]);
          nodes->SetSolution_Vel(iPoint_Local, iVar, Sol[iVar+nVar]);
          nodes->SetSolution_Vel_time_n(iPoint_Local, iVar, Sol[iVar+nVar]);
          nodes->SetSolution_Accel(iPoint_Local, iVar, Sol[iVar+2*nVar]);
          nodes->SetSolution_Accel_time_n(iPoint_Local, iVar, Sol[iVar+2*nVar]);
        }
        if (fluid_structure && !dynamic) {
          nodes->SetSolution_Pred(iPoint_Local, iVar, Sol[iVar]);
          nodes->SetSolution_Pred_Old(iPoint_Local, iVar, Sol[iVar]);
        }
        if (fluid_structure && discrete_adjoint){
          nodes->SetSolution_Old(iPoint_Local, iVar, Sol[iVar]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file. ---*/

  if (counter != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI. If dynamic, we also need to communicate the old solution. ---*/

  solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_FEA);
  solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_FEA);

  if (dynamic) {
    solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_FEA_OLD);
    solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_FEA_OLD);
  }
  if (fluid_structure && !dynamic) {
    solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_PRED);
    solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_PRED);

    solver[MESH_0][FEA_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION_PRED_OLD);
    solver[MESH_0][FEA_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION_PRED_OLD);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars; Restart_Vars = nullptr;
  delete [] Restart_Data; Restart_Data = nullptr;

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

void CFEASolver::FilterElementDensities(CGeometry *geometry, const CConfig *config)
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
  config->GetTopology_Optim_Projection(type,param);

  su2double *physical_rho = new su2double [nElement];

  /*--- "Rectify" the input, initialize the physical density with
  the design density (the filter function works in-place). ---*/
  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size))
  for (auto iElem=0ul; iElem<nElement; ++iElem) {
    su2double rho = element_properties[iElem]->GetDesignDensity();
    if      (rho > 1.0) physical_rho[iElem] = 1.0;
    else if (rho < 0.0) physical_rho[iElem] = 0.0;
    else                physical_rho[iElem] = rho;
  }

  geometry->FilterValuesAtElementCG(filter_radius, kernels, search_lim, physical_rho);

  SU2_OMP_PARALLEL
  {
    /*--- Apply projection. ---*/
    switch (type) {
      case NO_PROJECTION: break;
      case HEAVISIDE_UP:
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (auto iElem=0ul; iElem<nElement; ++iElem)
          physical_rho[iElem] = 1.0-exp(-param*physical_rho[iElem])+physical_rho[iElem]*exp(-param);
        break;
      case HEAVISIDE_DOWN:
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (auto iElem=0ul; iElem<nElement; ++iElem)
          physical_rho[iElem] = exp(-param*(1.0-physical_rho[iElem]))-(1.0-physical_rho[iElem])*exp(-param);
        break;
      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Unknown type of projection function",CURRENT_FUNCTION);
    }

    /*--- If input was out of bounds use the bound instead of the filtered
     value, useful to enforce solid or void regions (e.g. a skin). ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iElem=0ul; iElem<nElement; ++iElem) {
      su2double rho = element_properties[iElem]->GetDesignDensity();
      if      (rho > 1.0) element_properties[iElem]->SetPhysicalDensity(1.0);
      else if (rho < 0.0) element_properties[iElem]->SetPhysicalDensity(0.0);
      else element_properties[iElem]->SetPhysicalDensity(physical_rho[iElem]);
    }
  }

  delete [] physical_rho;

  /*--- For when this method is called directly, e.g. by the adjoint solver. ---*/
  topol_filter_applied = true;
}
