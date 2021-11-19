/*!
 * \file CSolverGradientSmoothing.cpp
 * \brief Main solver routines for the gradient smoothing problem.
 * \author T. Dick
 * \version 7.2.1 "Blackbird"
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

#include "../../include/solvers/CGradientSmoothingSolver.hpp"
#include "../../include/variables/CSobolevSmoothingVariable.hpp"
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

CGradientSmoothingSolver::CGradientSmoothingSolver(CGeometry *geometry, CConfig *config) : CSolver(false,true) {

  unsigned int iDim, jDim, marker_count=0;
  unsigned long iPoint;

  /*--- This solver does not perform recording operations --*/
  adjoint = false;

  /*--- general geometric settings ---*/
  nDim         = geometry->GetnDim();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nElement     = geometry->GetnElem();

  /*--- Here is where we assign the kind of each element ---*/

  /*--- First level: different possible terms of the equations  ---*/
  element_container = new CElement** [MAX_TERMS]();
  for (unsigned int iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS]();

  /*--- Initialize all subsequent levels ---*/
  for (unsigned int iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
    element_container[GRAD_TERM][iKind] = nullptr;
  }

  if (nDim == 2) {
    element_container[GRAD_TERM][EL_TRIA] = new CTRIA1();
    element_container[GRAD_TERM][EL_QUAD] = new CQUAD4();
    element_container[GRAD_TERM][EL_TRIA2] = new CTRIA3();
  }
  else if (nDim == 3) {
    element_container[GRAD_TERM][EL_TETRA] = new CTETRA1();
    element_container[GRAD_TERM][EL_HEXA]  = new CHEXA8();
    element_container[GRAD_TERM][EL_PYRAM] = new CPYRAM5();
    element_container[GRAD_TERM][EL_PRISM] = new CPRISM6();
    element_container[GRAD_TERM][EL_TETRA2] = new CTETRA4();
    element_container[GRAD_TERM][EL_PYRAM2] = new CPYRAM6();
  }

  /*--- for operations on surfaces we initalize the structures for nDim-1 ---*/
  if (config->GetSmoothOnSurface()) {
    if (nDim == 2) {
      element_container[GRAD_TERM][EL_LINE] = new CLINE();
    }
    else if (nDim == 3) {
      element_container[GRAD_TERM][EL_TRIA] = new CTRIA1();
      element_container[GRAD_TERM][EL_QUAD] = new CQUAD4();
      element_container[GRAD_TERM][EL_TRIA2] = new CTRIA3();
    }
  }

  Residual = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Residual[iDim] = 0.0;
  Solution = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Solution[iDim] = 0.0;
  mId_Aux    = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++){
    mId_Aux[iDim]    = new su2double[nDim];
  }
  for(iDim = 0; iDim < nDim; iDim++){
    for (jDim = 0; jDim < nDim; jDim++){
      mId_Aux[iDim][jDim]    = 0.0;
    }
    mId_Aux[iDim][iDim] = 1.0;
  }

  /*--- initializations for linear equation systems ---*/
  if ( !config->GetSmoothOnSurface() ) {
    if ( config->GetSepDim() ) {
      LinSysSol.Initialize(nPoint, nPointDomain, 1, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, 1, 0.0);
      Jacobian.Initialize(nPoint, nPointDomain, 1, 1, false, geometry, config, false, true);
    } else {
      LinSysSol.Initialize(nPoint, nPointDomain, nDim, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, nDim, 0.0);
      Jacobian.Initialize(nPoint, nPointDomain, nDim, nDim, false, geometry, config, false, true);
    }

  } else {
    if (config->GetSobMode() == PARAM_LEVEL_COMPLETE) {
      Jacobian.Initialize(nPoint, nPointDomain, nDim, nDim, false, geometry, config, false , true);
    } else {
      LinSysSol.Initialize(nPoint, nPointDomain, 1, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, 1, 0.0);
      Jacobian.Initialize(nPoint, nPointDomain, 1, 1, false, geometry, config, false, true);
    }
  }

  /*--- initializations for debugging ---*/
  if (config->GetSobMode()==DEBUG) {
    auxVec.Initialize(nPoint, nPointDomain, nDim, 1.0);
  }

  /*--- vectors needed for projection when working on complete system ---*/
  activeCoord.Initialize(nPoint, nPointDomain, nDim, 0.0);
  helperVecIn.Initialize(nPoint, nPointDomain, nDim, 0.0);
  helperVecOut.Initialize(nPoint, nPointDomain, nDim, 0.0);
  helperVecAux.Initialize(nPoint, nPointDomain, nDim, 0.0);

  /*--- Initialize the CVariable structure holding solution data ---*/
  nodes = new CSobolevSmoothingVariable(nPoint, nDim,  config);
  SetBaseClassPointerToNodes();

  /*--- Initialize the boundary of the boundary ---*/
  if (config->GetSmoothOnSurface()) {

    /*--- check which points are in more than one physical boundary ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (unsigned int iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
          long iVertex = geometry->nodes->GetVertex(iPoint, iMarker);
          if (iVertex >= 0) {
            marker_count++;
          }
        }
      }
      if (marker_count>=2) {
        nodes->MarkAsBoundaryPoint(iPoint);
      }
      marker_count = 0;
    }
  }

  /*--- Term ij of the Jacobian ---*/
  Jacobian_block = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Jacobian_block[iDim] = new su2double [nDim];
    for (jDim = 0; jDim < nDim; jDim++) {
      Jacobian_block[iDim][jDim] = 0.0;
    }
  }

  /*--- vector for the parameter gradient ---*/
  for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
    deltaP.push_back(0.0);
  }

  /*--- set the solver name for output purposes ---*/
  SolverName = "SOBOLEV";
}


CGradientSmoothingSolver::~CGradientSmoothingSolver(void) {

  unsigned int iDim;

  if (element_container != nullptr) {
    for (unsigned int iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (unsigned int jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }

  if (Jacobian_block != nullptr) {
    for (iDim = 0; iDim < nDim; ++iDim) {
      delete [] Jacobian_block[iDim];
    }
    delete [] Jacobian_block;
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    if (mId_Aux[iDim] != nullptr) delete [] mId_Aux[iDim];
  }
  if (mId_Aux != nullptr) delete [] mId_Aux;

  delete nodes;

}


void CGradientSmoothingSolver::ApplyGradientSmoothingVolume(CGeometry *geometry, CSolver *solver, CNumerics **numerics, CConfig *config) {

  /*--- current dimension if we run consecutive on each dimension ---*/
  dir = 0;

  /*--- Set vector and sparse matrix to 0 ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();

  Compute_StiffMatrix(geometry, numerics, config);

  if ( config->GetSepDim() ) {

    for (dir = 0; dir < nDim ; dir++) {

      Compute_Residual(geometry, solver, config);

      Impose_BC(geometry, numerics, config);

      Solve_Linear_System(geometry, config);

      WriteSensitivity(geometry, solver, config);

      LinSysSol.SetValZero();
      LinSysRes.SetValZero();
    }

  } else {

    Compute_Residual(geometry, solver, config);

    Impose_BC(geometry, numerics, config);

    Solve_Linear_System(geometry, config);

    WriteSensitivity(geometry, solver, config);

  }

}


void CGradientSmoothingSolver::ApplyGradientSmoothingSurface(CGeometry *geometry, CSolver *solver, CNumerics **numerics, CConfig *config, unsigned long val_marker) {

  /*--- Set vector and sparse matrix to 0 ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();

  Compute_Surface_StiffMatrix(geometry, numerics, config, val_marker);

  Compute_Surface_Residual(geometry, solver, config, val_marker);

  if ( config->GetDirichletSurfaceBound() ) {
    BC_Surface_Dirichlet(geometry, config, val_marker);
  }

  Solve_Linear_System(geometry, config);

  WriteSensitivity(geometry, solver, config, val_marker);

}


void CGradientSmoothingSolver::ApplyGradientSmoothingDV(CGeometry *geometry, CSolver *solver, CNumerics **numerics, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement, CConfig *config, vector<su2double> additionalGrad) {

  /// Set to 0
  Jacobian.SetValZero();

  /// record the parameterization
  if (rank == MASTER_NODE)  cout << " calculate the original gradient" << endl;
  RecordParameterizationJacobian(geometry, surface_movement, activeCoord, config);

  /// calculate the original gradinet
  CalculateOriginalGradient(geometry, grid_movement, config);

  /// if there was an initialization, add it too the original gradient
  if (additionalGrad.size()==deltaP.size()) {
    for (auto iDV = 0; iDV < deltaP.size(); iDV++) {
      deltaP[iDV] += additionalGrad[iDV];
    }
  }

  /// compute the Hessian column by column
  if (rank == MASTER_NODE)  cout << " computing the system matrix line by line" << endl;

  /// variable declarations
  unsigned nDVtotal=config->GetnDV_Total();
  unsigned column, row;
  vector<su2double> seedvector(nDVtotal, 0.0);
  hessian = MatrixType::Zero(nDVtotal, nDVtotal);

  /// get matrix vector product for this
  auto mat_vec = GetStiffnessMatrixVectorProduct(geometry, numerics, config);

  for (column=0; column<nDVtotal; column++) {

    if (rank == MASTER_NODE)  cout << "    working in column " << column << endl;

    /// create seeding
    std::fill(seedvector.begin(), seedvector.end(), 0.0);
    seedvector[column] = 1.0;

    // necessary for surface case!
    helperVecIn.SetValZero();
    helperVecOut.SetValZero();
    helperVecAux.SetValZero();

    /// forward projection
    ProjectDVtoMesh(geometry, seedvector, helperVecIn, activeCoord, config);

    /// matrix vector product in the middle
    if (config->GetSmoothOnSurface()) {

      /// perform MPI communication
      CSysMatrixComms::Initiate(helperVecIn, geometry, config, SOLUTION_MATRIX);
      CSysMatrixComms::Complete(helperVecIn, geometry, config, SOLUTION_MATRIX);

      /// compute the matrix vector product
      /// MatrixVector product operator does mpi comm at the end
      mat_vec(helperVecIn, helperVecAux);

      /// delete entries from halo cells, otherwise we get errors for parallel computations
      for (auto iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        for (auto iDim = 0; iDim < nDim; iDim++) {
          helperVecOut(iPoint,iDim) = helperVecAux(iPoint,iDim);
        }
      }

    } else {

      /// meshdeformation forward
      WriteVector2Geometry(geometry,config, helperVecIn);
      grid_movement->SetVolume_Deformation(geometry, config, false, true, true);
      ReadVector2Geometry(geometry,config, helperVecIn);

      /// perform MPI communication
      CSysMatrixComms::Initiate(helperVecIn, geometry, config, SOLUTION_MATRIX);
      CSysMatrixComms::Complete(helperVecIn, geometry, config, SOLUTION_MATRIX);

      /// compute the matrix vector product
      /// MatrixVector product operator does mpi comm at the end
      mat_vec(helperVecIn, helperVecAux);

      /// delete entries from halo cells, otherwise we get errors for parallel computations
      for (auto iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        for (auto iDim = 0; iDim < nDim; iDim++) {
          helperVecOut(iPoint,iDim) = helperVecAux(iPoint,iDim);
        }
      }

      /// mesh deformation backward
      WriteVector2Geometry(geometry,config, helperVecOut);
      grid_movement->SetVolume_Deformation(geometry, config, false, true, false);
      ReadVector2Geometry(geometry,config, helperVecOut);

    }

    /// reverse projection
    ProjectMeshToDV(geometry, helperVecOut, seedvector, activeCoord, config);

    /// extract projected direction
    hessian.col(column) = Eigen::Map<VectorType, Eigen::Unaligned>(seedvector.data(), seedvector.size());
  }

  /// output the matrix
  if (rank == MASTER_NODE) {
    ofstream SysMatrix(config->GetObjFunc_Hess_FileName());
    SysMatrix << hessian.format(CSVFormat);
    SysMatrix.close();
  }

  /// calculate and output the treated gradient
  QRdecomposition QR(hessian);
  VectorType b = Eigen::Map<VectorType, Eigen::Unaligned>(deltaP.data(), deltaP.size());
  VectorType x = QR.solve(b);
  deltaP = std::vector<su2double>(x.data(), x.data() + x.size());

  OutputDVGradient();

}


void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config){

  unsigned long iElem, iNode;
  unsigned int iDim, nNodes = 0, NelNodes, jNode;
  std::vector<unsigned long> indexNode(MAXNNODE_3D, 0.0);
  su2double val_Coord;
  int EL_KIND = 0;

  su2activematrix DHiDHj;
  su2double HiHj = 0.0;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element_container[GRAD_TERM][EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
      }

    }

    /*--- compute the contributions of the single elements inside the numerics container ---*/

    numerics[GRAD_TERM]->Compute_Tangent_Matrix(element_container[GRAD_TERM][EL_KIND], config);

    NelNodes = element_container[GRAD_TERM][EL_KIND]->GetnNodes();

    /*--- for all nodes add the contribution to the system Jacobian ---*/

    for (iNode = 0; iNode < NelNodes; iNode++) {

      for (jNode = 0; jNode < NelNodes; jNode++) {

        DHiDHj = element_container[GRAD_TERM][EL_KIND]->Get_DHiDHj(iNode, jNode);
        HiHj = element_container[GRAD_TERM][EL_KIND]->Get_HiHj(iNode, jNode);

        if ( config->GetSepDim() ) {

          Jacobian_block[0][0] = DHiDHj[dir][dir] + HiHj;
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_block);

        } else {

          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_block[iDim][iDim] = DHiDHj[iDim][iDim] + HiHj;
          }
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_block);

        }
      }
    }
  }
}


void CGradientSmoothingSolver::Compute_Surface_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config, unsigned long val_marker, unsigned int nSurfDim){

  unsigned long iElem, iPoint, iVertex, iSurfDim;
  unsigned int iNode, jNode, nNodes = 0, NelNodes;
  std::vector<unsigned long> indexNode(MAXNNODE_2D, 0.0);
  std::vector<unsigned long> indexVertex(MAXNNODE_2D, 0.0);
  int EL_KIND = 0;

  bool* visited = new bool[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
    visited[iPoint] = false;
  }

  su2activematrix DHiDHj;
  su2double HiHj = 0.0;

  su2activematrix Coord;

  /*--- Check if the current MPI rank has a part of the marker ---*/
  if (val_marker!=NOT_AVAILABLE) {

    /*--- Loops over all the elements ---*/
    for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

      /*--- Identify the kind of boundary element ---*/
      GetElemKindAndNumNodes(geometry->bound[val_marker][iElem]->GetVTK_Type(), EL_KIND, nNodes);

      /*--- Retrieve the boundary reference and current coordinates ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);
      }

      Coord = GetElementCoordinates(geometry, indexNode, EL_KIND);

      /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        for (iNode = 0; iNode < nNodes; iNode++) {
          if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
        }
      }

      /*--- compute the contributions of the single elements inside the numerics container ---*/
      numerics[GRAD_TERM]->SetCoord(Coord);
      numerics[GRAD_TERM]->Compute_Tangent_Matrix(element_container[GRAD_TERM][EL_KIND], config);

      NelNodes = element_container[GRAD_TERM][EL_KIND]->GetnNodes();

      /*--- for all nodes add the contribution to the system Jacobian ---*/

      for (iNode = 0; iNode < NelNodes; iNode++) {
        for (jNode = 0; jNode < NelNodes; jNode++) {

          DHiDHj = element_container[GRAD_TERM][EL_KIND]->Get_DHiDHj(iNode, jNode);
          HiHj = element_container[GRAD_TERM][EL_KIND]->Get_HiHj(iNode, jNode);

          /*--- TODO: need to do this for all block dimensions which is 2/3D respectively with the same component from surface 1/2D dimension! Maybe rework this inefficiency ---*/
          for (iSurfDim=0; iSurfDim<nSurfDim; iSurfDim++) {
            Jacobian_block[iSurfDim][iSurfDim] = DHiDHj[0][0] + HiHj;
          }
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_block);

          // store if we set something for a node
          if (visited[indexNode[iNode]]==false) { visited[indexNode[iNode]]=true; }
          if (visited[indexNode[jNode]]==false) { visited[indexNode[jNode]]=true; }

        }
      }
    }
  }

  // set the matrix to identity for non involved nodes
  // TODO: check if this should only be done for nodes in the domain?
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++){
    if (visited[iPoint]==false) {
      Jacobian.AddBlock(iPoint, iPoint, mId_Aux);
    }
  }
}


void CGradientSmoothingSolver::Compute_Residual(CGeometry *geometry, CSolver *solver, CConfig *config){

  unsigned long iElem;
  unsigned int iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  std::vector<unsigned long> indexNode(MAXNNODE_3D, 0.0);
  su2double Weight, Jac_X;

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        auto val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element_container[GRAD_TERM][EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
      }

    }

    element_container[GRAD_TERM][EL_KIND]->ClearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
    element_container[GRAD_TERM][EL_KIND]->ComputeGrad_Linear();
    unsigned int nGauss = element_container[GRAD_TERM][EL_KIND]->GetnGaussPoints();

    for (unsigned int iGauss = 0; iGauss < nGauss; iGauss++) {

      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      }

      Weight = element_container[GRAD_TERM][EL_KIND]->GetWeight(iGauss);
      Jac_X = element_container[GRAD_TERM][EL_KIND]->GetJ_X(iGauss);

      for (unsigned int iNode = 0; iNode < nNodes; iNode++) {

        if ( config->GetSepDim() ) {

          if (config->GetSobMode()==DEBUG) {
            Residual[dir] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * (auxVec.GetBlock(indexNode[iNode]))[dir];
            LinSysRes.AddBlock(indexNode[iNode], &Residual[dir]);
          } else {
            Residual[dir] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * nodes->GetSensitivity(indexNode[iNode], dir);
            LinSysRes.AddBlock(indexNode[iNode], &Residual[dir]);
          }

        } else {

          for (iDim = 0; iDim < nDim; iDim++) {

            if (config->GetSobMode()==DEBUG) {
              Residual[iDim] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * (auxVec.GetBlock(indexNode[iNode]))[iDim];
            } else {
              Residual[iDim] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * nodes->GetSensitivity(indexNode[iNode], iDim);
            }
          }
          LinSysRes.AddBlock(indexNode[iNode], Residual);

        }

        for (iDim = 0; iDim < nDim; iDim++) {
          Residual[iDim] = 0;
        }

      }
    }
  }
}


void CGradientSmoothingSolver::Compute_Surface_Residual(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned long val_marker){

  unsigned long iElem, iPoint, iVertex;
  unsigned int iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  std::vector<unsigned long> indexNode(MAXNNODE_2D, 0.0);
  std::vector<unsigned long> indexVertex(MAXNNODE_2D, 0.0);
  su2double Weight, Jac_X, normalSens = 0.0, norm;
  su2double* normal = NULL;
  su2activematrix Coord;

  /*--- Check if the current MPI rank has a part of the marker ---*/
  if (val_marker!=NOT_AVAILABLE) {

    for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

      /*--- Identify the kind of boundary element ---*/
      GetElemKindAndNumNodes(geometry->bound[val_marker][iElem]->GetVTK_Type(), EL_KIND, nNodes);

      /*--- Retrieve the boundary reference and current coordinates ---*/
      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);
      }

      Coord = GetElementCoordinates(geometry, indexNode, EL_KIND);

      /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        for (iNode = 0; iNode < nNodes; iNode++) {
          if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
        }
      }

      element_container[GRAD_TERM][EL_KIND]->ClearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
      element_container[GRAD_TERM][EL_KIND]->ComputeGrad_Linear(Coord);
      unsigned int nGauss = element_container[GRAD_TERM][EL_KIND]->GetnGaussPoints();

      for (unsigned int iGauss = 0; iGauss < nGauss; iGauss++) {

        Weight = element_container[GRAD_TERM][EL_KIND]->GetWeight(iGauss);
        Jac_X = element_container[GRAD_TERM][EL_KIND]->GetJ_X(iGauss);

        for (unsigned int iNode = 0; iNode < nNodes; iNode++) {

          normal = geometry->vertex[val_marker][indexVertex[iNode]]->GetNormal();
          norm = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            norm += normal[iDim]*normal[iDim];
          }
          norm = sqrt(norm);
          for (iDim = 0; iDim < nDim; iDim++) {
            normal[iDim] = normal[iDim] / norm;
          }

          for (iDim = 0; iDim < nDim; iDim++) {
            if (config->GetSobMode()==DEBUG) {
              normalSens += normal[iDim] * (auxVec.GetBlock(indexNode[iNode]))[iDim];
            } else {
              normalSens += normal[iDim] * nodes->GetSensitivity(indexNode[iNode], iDim);
            }
          }

          Residual[0] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * normalSens;
          LinSysRes.AddBlock(indexNode[iNode], Residual);

          Residual[0] = 0;
          normalSens = 0;

        }
      }
    }
  }
}


void CGradientSmoothingSolver::Impose_BC(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  unsigned int iMarker;

  /*--- Get the boundary markers and iterate over them ---------------------------------*/
  /* Notice that for no marker we automatically impose Zero Neumann boundary conditions */

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_SobolevBC(iMarker) == YES) {
      BC_Dirichlet(geometry, NULL, numerics, config, iMarker);
    }
  }

}


void CGradientSmoothingSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, unsigned int val_marker) {

  unsigned long iPoint, iVertex;
  const su2double zeros[MAXNDIM] = {0.0};

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Enforce solution, independed from dimension since zeros vector is long enough. ---*/
    LinSysSol.SetBlock(iPoint, zeros);
    LinSysRes.SetBlock(iPoint, zeros);
    Jacobian.EnforceSolutionAtNode(iPoint, zeros, LinSysRes);

  }

}


void CGradientSmoothingSolver::BC_Surface_Dirichlet(CGeometry *geometry, CConfig *config, unsigned int val_marker) {

  unsigned long iPoint, iVertex;
  const su2double zeros[MAXNDIM] = {0.0};

  /*--- Check if the current MPI rank has a part of the marker ---*/
  if (val_marker!=NOT_AVAILABLE) {
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

      /*--- Get node index ---*/
      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if ( nodes->IsBoundaryPoint(iPoint) ) {

        LinSysSol.SetBlock(iPoint, zeros);
        LinSysRes.SetBlock(iPoint, zeros);
        Jacobian.EnforceSolutionAtNode(iPoint, zeros, LinSysRes);

      }
    }
  }
}


void CGradientSmoothingSolver::Solve_Linear_System(CGeometry *geometry, CConfig *config){

  unsigned long IterLinSol = 0;

  /* For MPI prescribe vector entries across the ranks before solving the system.
   * Analog to FEA solver this is only done for the solution */
  CSysMatrixComms::Initiate(LinSysSol, geometry, config);
  CSysMatrixComms::Complete(LinSysSol, geometry, config);

  // elimination shedule
  Set_VertexEliminationSchedule(geometry, config);

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  SU2_OMP_MASTER
  {
    SetIterLinSolver(IterLinSol);
    SetResLinSolver(System.GetResidual());
  }

}


CSysMatrixVectorProduct<su2matvecscalar> CGradientSmoothingSolver::GetStiffnessMatrixVectorProduct(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  bool surf = config->GetSmoothOnSurface();

  /*--- Compute the sparse stiffness matrix ---*/
  if (surf) {
    unsigned long dvMarker=NOT_AVAILABLE;
    for (unsigned int iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ( config->GetMarker_All_DV(iMarker) == YES ) {
        dvMarker=iMarker;
      }
    }
    Compute_Surface_StiffMatrix(geometry, numerics, config, dvMarker, nDim);
  } else {
    Compute_StiffMatrix(geometry, numerics, config);
  }

  return CSysMatrixVectorProduct<su2matvecscalar>(Jacobian, geometry, config);
}


void CGradientSmoothingSolver::CalculateOriginalGradient(CGeometry *geometry, CVolumetricMovement *grid_movement, CConfig *config) {

  if (rank == MASTER_NODE) cout << endl << "Calculating the original DV gradient." << endl;

  WriteSens2Geometry(geometry,config);

  grid_movement->SetVolume_Deformation(geometry, config, false, true);

  ReadVector2Geometry(geometry,config, helperVecOut);

  ProjectMeshToDV(geometry, helperVecOut, deltaP, activeCoord, config);

  OutputDVGradient("orig_grad.dat");
}


void CGradientSmoothingSolver::RecordTapeAndCalculateOriginalGradient(CGeometry *geometry, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement, CConfig *config) {

  /// record the parameterization
  if (rank == MASTER_NODE)  cout << " calculate the original gradient" << endl;
  RecordParameterizationJacobian(geometry, surface_movement, activeCoord, config);

  /// calculate the original gradinet
  CalculateOriginalGradient(geometry, grid_movement, config);

}


void CGradientSmoothingSolver::OutputDVGradient(string out_file) {

  unsigned iDV;
  if (rank == MASTER_NODE) {
    ofstream delta_p (out_file);
    delta_p.precision(17);
    for (iDV = 0; iDV < deltaP.size(); iDV++) {
      delta_p << deltaP[iDV] << ",";
    }
    delta_p.close();
  }

}


void CGradientSmoothingSolver::RecordParameterizationJacobian(CGeometry *geometry, CSurfaceMovement *surface_movement, CSysVector<su2double>& registeredCoord, CConfig *config) {

  unsigned int nDim, nMarker, nDV, nDV_Value, nPoint, nVertex;
  unsigned int iDV, iDV_Value, iMarker, iPoint, iVertex, iDim;
  su2double* VarCoord;
  su2double** DV_Value;

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nPoint  = geometry->GetnPoint();
  nDV     = config->GetnDV();

  /*--- Start recording of operations ---*/

  AD::Reset();

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  DV_Value = config->GetDV_Pointer();
  for (iDV = 0; iDV < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization of su2double with 0.0 resets the index ---*/
      DV_Value[iDV][iDV_Value] = 0.0;
      AD::RegisterInput(DV_Value[iDV][iDV_Value]);

    }
  }

  /*--- Call the surface deformation routine ---*/
  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Register Outputs --- */
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim=0; iDim<nDim; iDim++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          registeredCoord(iPoint,iDim) = VarCoord[iDim];
        }
      }
    }
  }
  for (iPoint = 0; iPoint<nPoint; iPoint++) {
    for (iDim=0; iDim<nDim; iDim++) {
      AD::RegisterOutput(registeredCoord(iPoint,iDim));
    }
  }

  /*--- Stop the recording --- */
  AD::StopRecording();

}


void CGradientSmoothingSolver::ProjectDVtoMesh(CGeometry *geometry, std::vector<su2double>& seeding, CSysVector<su2matvecscalar>& result, CSysVector<su2double>& registeredCoord, CConfig *config) {

  unsigned int nDim, nMarker, nDV, nDV_Value, nVertex;
  unsigned int iDV, iDV_Value, iDV_index, iMarker, iVertex, iPoint, iDim;
  su2double** DV_Value = config->GetDV_Pointer();

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  /*--- Seeding for the DV. --*/
  iDV_index = 0;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      SU2_TYPE::SetDerivative(DV_Value[iDV][iDV_Value], SU2_TYPE::GetValue(seeding[iDV_index]));
      iDV_index++;
    }
  }

  AD::ComputeAdjointForward();

  /*--- Extract sensitivities ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim=0; iDim<nDim; iDim++) {
          result(iPoint,iDim) = SU2_TYPE::GetDerivative(registeredCoord(iPoint,iDim));
        }
      }
    }
  }

  AD::ClearAdjoints();

}


void CGradientSmoothingSolver::ProjectMeshToDV(CGeometry *geometry, CSysVector<su2matvecscalar>& sensitivity, std::vector<su2double>& output, CSysVector<su2double>& registeredCoord, CConfig *config) {

  /*--- adjoint surface deformation ---*/

  unsigned int nDim, nMarker, nDV, nDV_Value, nVertex;
  unsigned int iDV, iDV_Value, iDV_index, iPoint, iDim, iMarker, iVertex;
  su2double** DV_Value = config->GetDV_Pointer();
  su2double my_Gradient, localGradient;

  // get some numbers from config
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  /*--- Set the seeding appropriately ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++){
          SU2_TYPE::SetDerivative(registeredCoord(iPoint,iDim), SU2_TYPE::GetValue(sensitivity(iPoint,iDim)));
        }
      }
    }
  }

  /*--- Compute derivatives and extract gradient ---*/
  AD::ComputeAdjoint();

  iDV_index = 0;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      my_Gradient = SU2_TYPE::GetDerivative(DV_Value[iDV][iDV_Value]);
      #ifdef HAVE_MPI
        SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #else
        localGradient = my_Gradient;
      #endif
      output[iDV_index] = localGradient;
      iDV_index++;
    }
  }

  AD::ClearAdjoints();

}


void CGradientSmoothingSolver::SetSensitivity(CGeometry *geometry, CConfig *config, CSolver *solver) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      nodes->SetSensitivity(iPoint,iDim, solver->GetNodes()->GetSensitivity(iPoint,iDim));
    }
  }
}


void CGradientSmoothingSolver::OutputSensitivity(CGeometry *geometry, CConfig *config, CSolver *solver) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      solver->GetNodes()->SetSensitivity(iPoint,iDim, nodes->GetSensitivity(iPoint,iDim));
    }
  }
}


void CGradientSmoothingSolver::WriteSensitivity(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned long val_marker){

  unsigned long iPoint, total_index;
  unsigned int iDim;
  su2double* normal;
  su2double norm;

  /*--- split between surface and volume first, to avoid mpi ranks with no part of the marker to write back nonphysical solutions in the surface case ---*/
  if ( config->GetSmoothOnSurface() ) {
    if( val_marker!=NOT_AVAILABLE ) {
      for (unsigned long iVertex =0; iVertex<geometry->nVertex[val_marker]; iVertex++)  {
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        normal = geometry->vertex[val_marker][iVertex]->GetNormal();
        norm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          norm += normal[iDim]*normal[iDim];
        }
        norm = sqrt(norm);
        for (iDim = 0; iDim < nDim; iDim++) {
          normal[iDim] = normal[iDim] / norm;
        }
        for (iDim = 0; iDim < nDim; iDim++) {
          this->GetNodes()->SetSensitivity(iPoint, iDim, normal[iDim]*LinSysSol[iPoint]);
        }
      }
    }
  } else {
    if ( config->GetSepDim() ) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        this->GetNodes()->SetSensitivity(iPoint, dir, LinSysSol[iPoint]);
      }
    } else {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint*nDim + iDim;
          this->GetNodes()->SetSensitivity(iPoint, iDim, LinSysSol[total_index]);
        }
      }
    }
  }
}


void CGradientSmoothingSolver::WriteSens2Geometry(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->SetSensitivity(iPoint,iDim, nodes->GetSensitivity(iPoint,iDim));
    }
  }
}


void CGradientSmoothingSolver::ReadSens2Geometry(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      nodes->SetSensitivity(iPoint, iDim, geometry->GetSensitivity(iPoint,iDim));
    }
  }
}


void CGradientSmoothingSolver::WriteVector2Geometry(CGeometry *geometry, CConfig *config, CSysVector<su2matvecscalar>& vector) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->SetSensitivity(iPoint,iDim, vector[iPoint*nDim+iDim]);
    }
  }
}


void CGradientSmoothingSolver::ReadVector2Geometry(CGeometry *geometry, CConfig *config, CSysVector<su2matvecscalar> &vector) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      vector[iPoint*nDim+iDim] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint,iDim));
    }
  }
}


su2activematrix CGradientSmoothingSolver::GetElementCoordinates(CGeometry *geometry, std::vector<unsigned long>& indexNode, int EL_KIND) {

  su2activematrix Coord;

  switch (EL_KIND) {

  case EL_LINE:

    Coord.resize(2,2);
    for(auto iNode=0; iNode<2; iNode++) {
      for(auto iDim=0; iDim<2; iDim++) {
        Coord[iNode][iDim] = Get_ValCoord(geometry, indexNode[iNode], iDim);
      }
    }
    break;

  case EL_TRIA:

    Coord.resize(3,3);
    for(auto iNode=0; iNode<3; iNode++) {
      for(auto iDim=0; iDim<3; iDim++) {
        Coord[iNode][iDim] = Get_ValCoord(geometry, indexNode[iNode], iDim);
      }
    }
    break;

  case EL_TRIA2:

    Coord.resize(3,3);
    for(auto iNode=0; iNode<3; iNode++) {
      for(auto iDim=0; iDim<3; iDim++) {
        Coord[iNode][iDim] = Get_ValCoord(geometry, indexNode[iNode], iDim);
      }
    }
    break;

  case EL_QUAD:

    Coord.resize(4,3);
    for(auto iNode=0; iNode<4; iNode++) {
      for(auto iDim=0; iDim<3; iDim++) {
        Coord[iNode][iDim] = Get_ValCoord(geometry, indexNode[iNode], iDim);
      }
    }
    break;

  default:
    std::cout << "Type of element is not supported. " <<std::endl;

  }

  return Coord;

}


void CGradientSmoothingSolver::Set_VertexEliminationSchedule(CGeometry *geometry, CConfig *config) {

  /*--- Store global point indices of essential BC markers. ---*/
  vector<unsigned long> myPoints;
  unsigned long iPoint;

  /*--- get global indexes of Dirichlet boundaries ---*/

  if (config->GetSmoothOnSurface()) {

    /*--- Find Marker_DV if this node has any ---*/
    unsigned long dvMarker=NOT_AVAILABLE;
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES) {
        dvMarker = iMarker;
      }
    }

    /*--- surface case: design boundary border if Dirichlet condition is set and rank holds part of it. ---*/
    if ( config->GetDirichletSurfaceBound() && dvMarker!=NOT_AVAILABLE) {
      for (auto iVertex = 0; iVertex < geometry->nVertex[dvMarker]; iVertex++) {
        /*--- Get node index ---*/
        iPoint = geometry->vertex[dvMarker][iVertex]->GetNode();
        if ( nodes->IsBoundaryPoint(iPoint) ) {
          myPoints.push_back(geometry->nodes->GetGlobalIndex(iPoint));
        }
      }
    }

  } else {

    /*--- volume case: all boundaries where Dirichlet condition is set. ---*/
    vector<unsigned short> markers;
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_SobolevBC(iMarker) == YES) {
        markers.push_back(iMarker);
      }
    }
    for (auto iMarker : markers) {
      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        myPoints.push_back(geometry->nodes->GetGlobalIndex(iPoint));
      }
    }
  }

  /*--- communicate the boundary points ---*/

  const unordered_set<unsigned long> markerPoints(myPoints.begin(), myPoints.end());

  vector<unsigned long> numPoints(size);
  unsigned long num = myPoints.size();
  SU2_MPI::Allgather(&num, 1, MPI_UNSIGNED_LONG, numPoints.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Global to local map for the halo points of the rank (not covered by the CGeometry map). ---*/
  unordered_map<unsigned long, unsigned long> Global2Local;
  for (auto iPoint = nPointDomain; iPoint < nPoint; ++iPoint) {
    Global2Local[geometry->nodes->GetGlobalIndex(iPoint)] = iPoint;
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

  /*--- eliminate the extra vertices ---*/

  for (auto iPoint : ExtraVerticesToEliminate) {
    Jacobian.EnforceSolutionAtNode(iPoint, LinSysSol.GetBlock(iPoint), LinSysRes);
  }

}
