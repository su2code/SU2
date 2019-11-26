/*!
 * \file solver_gradient_smoothing.cpp
 * \brief Main subroutines for the gradient smoothing problem.
 * \author T. Dick
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
#include "../../Common/include/adt_structure.hpp"
#include "../include/variables/CSobolevSmoothingVariable.hpp"
#include <algorithm>


CGradientSmoothingSolver::CGradientSmoothingSolver(CGeometry *geometry, CConfig *config) : CSolver(false,true) {

  unsigned short iDim, jDim;
  unsigned int marker_count=0;
  unsigned long iPoint;

  /*--- general geometric settings ---*/
  nDim         = geometry->GetnDim();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nElement     = geometry->GetnElem();

  /*--- Element container structure ---*/

  /*--- First level: only the FEA_TERM is considered ---*/
  element_container = new CElement** [12];
  element_container[GRAD_TERM] = new CElement* [MAX_FE_KINDS];

  /*--- Initialize all subsequent levels ---*/
  for (unsigned short iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
    element_container[GRAD_TERM][iKind] = NULL;
  }

  if (nDim == 2) {
    element_container[GRAD_TERM][EL_TRIA] = new CTRIA1(nDim, config);
    element_container[GRAD_TERM][EL_QUAD] = new CQUAD4(nDim, config);
    if (config->GetSecOrdQuad()) {
      element_container[GRAD_TERM][EL_TRIA2] = new CTRIA3(nDim, config);
    }
  }
  else if (nDim == 3) {
    element_container[GRAD_TERM][EL_TETRA] = new CTETRA1(nDim, config);
    element_container[GRAD_TERM][EL_HEXA]  = new CHEXA8(nDim, config);
    element_container[GRAD_TERM][EL_PYRAM] = new CPYRAM5(nDim, config);
    element_container[GRAD_TERM][EL_PRISM] = new CPRISM6(nDim, config);
    if (config->GetSecOrdQuad()) {
      element_container[GRAD_TERM][EL_TETRA2] = new CTETRA4(nDim, config);
      element_container[GRAD_TERM][EL_PYRAM2] = new CPYRAM6(nDim, config);
    }
  }

  /*--- for operations on surfaces we initalize the structures for nDim-1 ---*/
  if (config->GetSmoothOnSurface()) {
    if (nDim == 2) {
      element_container[GRAD_TERM][EL_LINE] = new CLINE(nDim-1, config);
    }
    else if (nDim == 3) {
      element_container[GRAD_TERM][EL_TRIA] = new CTRIA1(nDim-1, config);
      element_container[GRAD_TERM][EL_QUAD] = new CQUAD4(nDim-1, config);
      if (config->GetSecOrdQuad()) {
        element_container[GRAD_TERM][EL_TRIA2] = new CTRIA3(nDim-1, config);
      }
    }
  }

  Residual = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Residual[iDim] = 0.0;
  Solution = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Solution[iDim] = 0.0;
  mZeros_Aux = new su2double *[nDim];
  mId_Aux    = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++){
    mZeros_Aux[iDim] = new su2double[nDim];
    mId_Aux[iDim]    = new su2double[nDim];
  }
  for(iDim = 0; iDim < nDim; iDim++){
    for (jDim = 0; jDim < nDim; jDim++){
      mZeros_Aux[iDim][jDim] = 0.0;
      mId_Aux[iDim][jDim]    = 0.0;
    }
    mId_Aux[iDim][iDim] = 1.0;
  }

  /*--- linear system ---*/
  if ( !config->GetSmoothOnSurface() ) {
    if ( config->GetSepDim() ) {
      LinSysSol.Initialize(nPoint, nPointDomain, 1, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, 1, 0.0);
      Jacobian.Initialize(nPoint, nPointDomain, 1, 1, false, geometry, config);
    } else {
      LinSysSol.Initialize(nPoint, nPointDomain, nDim, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, nDim, 0.0);
      Jacobian.Initialize(nPoint, nPointDomain, nDim, nDim, false, geometry, config);
    }

    // initialize auxiliar helper vectors
    auxVecInp.Initialize(nPoint, nPointDomain, nDim, 0.0);
    auxVecRHS.Initialize(nPoint, nPointDomain, nDim, 0.0);
    auxVecOut.Initialize(nPoint, nPointDomain, nDim, 0.0);
  }

  /*--- Initialize the boundary of the boundary ---*/

  if (config->GetSmoothOnSurface()) {
    nodes = new CSobolevSmoothingVariable(nPoint, nDim,  config);
    SetBaseClassPointerToNodes();

    /*--- check which points are in more than one physical boundary ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {


      for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
        if (iVertex >= 0) {
          marker_count++;
        }
      }
      if (marker_count>=2) {
        nodes->MarkAsBoundaryPoint(iPoint);
      }
      marker_count = 0;
    }
  }

  /*--- Term ij of the Jacobian ---*/
  Jacobian_ij = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Jacobian_ij[iDim] = new su2double [nDim];
    for (jDim = 0; jDim < nDim; jDim++) {
      Jacobian_ij[iDim][jDim] = 0.0;
    }
  }


}


CGradientSmoothingSolver::~CGradientSmoothingSolver(void) {

  unsigned short iDim;

  if (element_container != NULL) {
    for (unsigned short iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (unsigned short jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        if (element_container[iVar][jVar] != NULL) delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }

  if (Residual != NULL) delete [] Residual;
  if (Solution != NULL) delete [] Solution;

  for (iDim = 0; iDim < nDim; iDim++) {
    if (mZeros_Aux[iDim] != NULL) delete [] mZeros_Aux[iDim];
    if (mId_Aux[iDim] != NULL) delete [] mId_Aux[iDim];
  }
  if (mZeros_Aux != NULL) delete [] mZeros_Aux;
  if (mId_Aux != NULL) delete [] mId_Aux;

  if (nodes != NULL) delete nodes;

}


void CGradientSmoothingSolver::ApplyGradientSmoothing(CGeometry *geometry, CSolver *solver, CNumerics **numerics, CConfig *config) {

  dir = 0;

  /*--- Initialize vector and sparse matrix ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();

  Compute_StiffMatrix(geometry, numerics, config);

  ofstream matrix ("matrix.dat");
  Jacobian.printMat(matrix);
  matrix.close();

  if ( config->GetSepDim() ) {

    for (dir = 0; dir < nDim ; dir++) {

      for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
        auxVecInp.SetBlock(iPoint, dir, solver->GetNodes()->GetSensitivity(iPoint ,dir));
      }

      Compute_Residual(geometry, solver, config);
      Impose_BC(geometry, numerics, config);

      for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
        auxVecRHS.SetBlock(iPoint, dir, LinSysRes.GetBlock(iPoint,0));
      }

      Solve_Linear_System(geometry, config);
      Set_Sensitivities(geometry, solver, config);

      ofstream result ("result.txt");
      LinSysSol.printVec(result);
      result.close();

      for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
        auxVecOut.SetBlock(iPoint, dir, solver->GetNodes()->GetSensitivity(iPoint ,dir));
      }

      LinSysSol.SetValZero();
      LinSysRes.SetValZero();
    }

    ofstream input ("input.txt");
    auxVecInp.printVec(input);
    input.close();

    ofstream rhs ("rhs.txt");
    auxVecRHS.printVec(rhs);
    rhs.close();

    ofstream output ("output.txt");
    auxVecOut.printVec(output);
    output.close();

  } else {

    for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
      su2double x = geometry->node[iPoint]->GetCoord(0);
      su2double y = geometry->node[iPoint]->GetCoord(1);
      for (auto iDim = 0; iDim < nDim ; iDim++) {
        auxVecInp.SetBlock(iPoint, iDim, -2.0*3.14159265359*3.14159265359*sin(3.14159265359*x)*sin(3.14159265359*y) );
      }
    }

    /*

    cos(3.14159265359*geometry->node[iPoint]->GetCoord(0))*cos(3.14159265359*geometry->node[iPoint]->GetCoord(1))

    */

    Compute_Residual(geometry, solver, config);

    Impose_BC(geometry, numerics, config);

    Solve_Linear_System(geometry, config);

    Set_Sensitivities(geometry, solver, config);

    for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
      for (auto iDim = 0; iDim < nDim ; iDim++) {
        auxVecOut.SetBlock(iPoint, iDim, solver->GetNodes()->GetSensitivity(iPoint ,iDim));
      }
    }

    ofstream input ("input.txt");
    auxVecInp.printVec(input);
    input.close();

    ofstream rhs ("rhs.txt");
    LinSysRes.printVec(rhs);
    rhs.close();

    ofstream output ("output.txt");
    auxVecOut.printVec(output);
    output.close();

  }

}


void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config){

  unsigned long iElem;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;

  su2double **DHiDHj = NULL;
  su2double HiHj = 0.0;

  unsigned short NelNodes, jNode;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_PYRAM;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_PRISM;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}

    // if we need higher order quadrature rules overide some of the element kinds
    if (config->GetSecOrdQuad()) {
      if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)    {nNodes = 3; EL_KIND = EL_TRIA2;}
      if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {nNodes = 4; EL_KIND = EL_TETRA2;}
      if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)     {nNodes = 6; EL_KIND = EL_PYRAM2;}
    }

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element_container[GRAD_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
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

          Jacobian_ij[0][0] = DHiDHj[dir][dir] + HiHj;
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

        } else {

          for (iDim = 0; iDim < nDim; iDim++) {
            DHiDHj[iDim][iDim] += HiHj;
          }
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], DHiDHj);

        }
      }
    }
  }
}


void CGradientSmoothingSolver::Compute_Surface_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config, unsigned long val_marker){

  // TO DO
  // - check initialisation of nElem_Bound and bound

  unsigned long iElem, iPoint, iVertex;
  unsigned short iNode, jNode, iDim, nNodes = 0, NelNodes;
  std::vector<unsigned long> indexNode(8, 0.0);
  std::vector<unsigned long> indexVertex(8, 0.0);
  su2double val_Coord;
  int EL_KIND = 0;

  su2double **DHiDHj = NULL;
  su2double HiHj = 0.0;

  std::vector<std::vector<su2double>> Coord;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    /*--- Identify the kind of boundary element ---*/
    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == LINE)           {nNodes = 2; EL_KIND = EL_LINE;}
    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL)  {nNodes = 4; EL_KIND = EL_QUAD;}

    if (config->GetSecOrdQuad()) {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       {nNodes = 3; EL_KIND = EL_TRIA2;}
    }

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

        DHiDHj[0][0] += HiHj;

        auto meshVertexI = geometry->node[indexNode[iNode]]->GetGlobalIndex();
        auto meshVertexJ = geometry->node[indexNode[jNode]]->GetGlobalIndex();
        std::cout << "Setting block " << indexVertex[iNode] << ", " << indexVertex[jNode] << " aka " << meshVertexI << ", " << meshVertexJ << " to: " << DHiDHj[0][0] << std::endl;
        Jacobian.AddBlock(indexVertex[iNode], indexVertex[jNode], DHiDHj);

      }
    }
  }
}


void CGradientSmoothingSolver::Compute_Residual(CGeometry *geometry, CSolver *solver, CConfig *config){

  unsigned long iElem;
  unsigned short iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  std::vector<unsigned long> indexNode(8, 0.0);
  su2double Weight, Jac_X;

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_PYRAM;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_PRISM;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}

    // if we need higher order quadrature rules overide some of the element kinds
    if (config->GetSecOrdQuad()) {
      if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)    {nNodes = 3; EL_KIND = EL_TRIA2;}
      if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {nNodes = 4; EL_KIND = EL_TETRA2;}
      if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)     {nNodes = 6; EL_KIND = EL_PYRAM2;}
    }

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        auto val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element_container[GRAD_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }

    }

    element_container[GRAD_TERM][EL_KIND]->clearElement(true);       /*--- Restarts the element: avoids adding over previous results in other elements --*/
    element_container[GRAD_TERM][EL_KIND]->ComputeGrad_Linear();
    unsigned short nGauss = element_container[GRAD_TERM][EL_KIND]->GetnGaussPoints();

    for (unsigned short iGauss = 0; iGauss < nGauss; iGauss++) {

      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      }

      Weight = element_container[GRAD_TERM][EL_KIND]->GetWeight(iGauss);
      Jac_X = element_container[GRAD_TERM][EL_KIND]->GetJ_X(iGauss);

      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {

        if ( config->GetSepDim() ) {

          if (config->GetSobDebugMode()) {
            Residual[dir] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * auxVecInp.GetBlock(indexNode[iNode], dir);
            LinSysRes.AddBlock(indexNode[iNode], &Residual[dir]);
          } else {
            Residual[dir] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * solver->GetNodes()->GetSensitivity(indexNode[iNode], dir);
            LinSysRes.AddBlock(indexNode[iNode], &Residual[dir]);
          }

        } else {

          for (iDim = 0; iDim < nDim; iDim++) {

            if (config->GetSobDebugMode()) {
              Residual[iDim] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * auxVecInp.GetBlock(indexNode[iNode], iDim);
            } else {
              Residual[iDim] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * solver->GetNodes()->GetSensitivity(indexNode[iNode], iDim);
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
  unsigned short iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  std::vector<unsigned long> indexNode(8, 0.0);
  std::vector<unsigned long> indexVertex(8, 0.0);
  su2double Weight, Jac_X, normalSens = 0.0, norm;
  su2double* normal = NULL;
  std::vector<std::vector<su2double>> Coord;

  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    /*--- Identify the kind of boundary element ---*/
    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == LINE)           {nNodes = 2; EL_KIND = EL_LINE;}
    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL)  {nNodes = 4; EL_KIND = EL_QUAD;}

    if (config->GetSecOrdQuad()) {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       {nNodes = 3; EL_KIND = EL_TRIA2;}
    }

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

    /*--- Retrieve the reference normal for one of the points. They go INSIDE the structural domain. ---*/
    normal = geometry->vertex[val_marker][indexVertex[0]]->GetNormal();
    norm = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      norm += normal[iDim]*normal[iDim];
    }
    norm = sqrt(norm);
    for (iDim = 0; iDim < nDim; iDim++) {
      normal[iDim] = normal[iDim] / norm;
    }

    element_container[GRAD_TERM][EL_KIND]->clearElement(true);       /*--- Restarts the element: avoids adding over previous results in other elements --*/
    element_container[GRAD_TERM][EL_KIND]->ComputeGrad_Linear(Coord);
    unsigned short nGauss = element_container[GRAD_TERM][EL_KIND]->GetnGaussPoints();

    for (unsigned short iGauss = 0; iGauss < nGauss; iGauss++) {

      Weight = element_container[GRAD_TERM][EL_KIND]->GetWeight(iGauss);
      Jac_X = element_container[GRAD_TERM][EL_KIND]->GetJ_X(iGauss);

      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {

        for (iDim = 0; iDim < nDim; iDim++) {
          if (config->GetSobDebugMode()) {
            normalSens += normal[iDim] * auxVecInp.GetBlock(indexVertex[iNode], iDim);
          } else {
            /*--- use indexNode here since we want the sensitivity (3D info) ---*/
            normalSens += normal[iDim] * solver->GetNodes()->GetSensitivity(indexNode[iNode], iDim);
          }
        }
        Residual[0] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * normalSens;
        LinSysRes.AddBlock(indexVertex[iNode], Residual);

        Residual[0] = 0;
        normalSens = 0;

      }
    }
  }
}



void CGradientSmoothingSolver::Impose_BC(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  unsigned short iMarker;
  /*--- Get the boundary markers and iterate over them ---------------------------------*/
  /* Notice that for no marker we automatically impose Zero Neumann boundary conditions */

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_SobolevBC(iMarker) == YES) {
      BC_Dirichlet(geometry, NULL, numerics, config, iMarker);
    }
  }

}


void CGradientSmoothingSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, unsigned short val_marker) {


  unsigned long iPoint, iVertex;
  unsigned long iVar, jVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if ( config->GetSepDim() ) {

      su2double one = 1.0;
      su2double zero = 0.0;

      if (geometry->node[iPoint]->GetDomain()) {

        LinSysRes.SetBlock(iPoint, &zero);
        LinSysSol.SetBlock(iPoint, &zero);

        for (iVar = 0; iVar < nPoint; iVar++) {
          if (iVar==iPoint) {
            Jacobian.SetBlock(iVar,iPoint, &one);
          }
          else {
            Jacobian.SetBlock(iVar,iPoint, &zero);
          }
        }
        /*--- Delete the rows for a particular node ---*/
        for (jVar = 0; jVar < nPoint; jVar++) {
          if (iPoint!=jVar) {
            Jacobian.SetBlock(iPoint,jVar, &zero);
          }
        }

      } else {
        /*--- Delete the column (iPoint is halo so Send/Recv does the rest) ---*/
        for (iVar = 0; iVar < nPoint; iVar++) {
          Jacobian.SetBlock(iVar,iPoint, &zero);
        }
      }

    } else {

      if (geometry->node[iPoint]->GetDomain()) {

       if (nDim == 2) {
          Solution[0] = 0.0;  Solution[1] = 0.0;
          Residual[0] = 0.0;  Residual[1] = 0.0;
        }
        else {
          Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
          Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
        }

        LinSysRes.SetBlock(iPoint, Residual);
        LinSysSol.SetBlock(iPoint, Solution);

        /*--- STRONG ENFORCEMENT OF THE DIRICHLET BOUNDARY CONDITION ---*/
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

      } else {
        /*--- Delete the column (iPoint is halo so Send/Recv does the rest) ---*/
        for (iVar = 0; iVar < nPoint; iVar++) Jacobian.SetBlock(iVar,iPoint,mZeros_Aux);
      }

    }

  }

}


void CGradientSmoothingSolver::BC_Surface_Dirichlet(CGeometry *geometry, CConfig *config, unsigned short val_marker) {


  unsigned long iPoint, iVertex;
  unsigned long iVar, jVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if ( nodes->IsBoundaryPoint(iPoint) ) {

      std::cout << "Node " << iPoint << " at ("<< geometry->node[iPoint]->GetCoord(0) << ", " << geometry->node[iPoint]->GetCoord(1) << ") is a bound of a bound." << std::endl;

      su2double one = 1.0;
      su2double zero = 0.0;

      if (geometry->node[iPoint]->GetDomain()) {

        LinSysRes.SetBlock(iVertex, &zero);
        LinSysSol.SetBlock(iVertex, &zero);

        for (iVar = 0; iVar < geometry->nVertex[val_marker]; iVar++) {
          if (iVar==iVertex) {
            std::cout << "Setting block " << iVar << ", " << iVertex << std::endl;
            Jacobian.SetBlock(iVar,iVertex, &one);
          }
          else {
            std::cout << "Setting block " << iVar << ", " << iVertex << std::endl;
            Jacobian.SetBlock(iVar,iVertex, &zero);
          }
        }
        /*--- Delete the rows for a particular node ---*/
        for (jVar = 0; jVar < geometry->nVertex[val_marker]; jVar++) {
          if (iVertex!=jVar) {
            std::cout << "Setting block " << iVertex << ", " << jVar << std::endl;
            Jacobian.SetBlock(iVertex,jVar, &zero);
          }
        }

      } else {
        /*--- Delete the column (iPoint is halo so Send/Recv does the rest) ---*/
        for (iVar = 0; iVar < geometry->nVertex[val_marker]; iVar++) {
          std::cout << "Setting block " << iVar<< ", " << iVertex << std::endl;
          Jacobian.SetBlock(iVar,iVertex, &zero);
        }
      }
    }
  }

}


// For now: left empty since there is no calculation necessary for zero Neumann boundaries.
void CGradientSmoothingSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, unsigned short val_marker) {

}


void CGradientSmoothingSolver::Solve_Linear_System(CGeometry *geometry, CConfig *config){

  unsigned long IterLinSol = 0;

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  SetIterLinSolver(IterLinSol);

}


void CGradientSmoothingSolver::Set_Sensitivities(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned long val_marker){

  unsigned long iPoint, total_index;
  unsigned short iDim;

  if ( config->GetSepDim() ) {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      solver->GetNodes()->SetSensitivity(iPoint, dir, LinSysSol[iPoint]);
    }

  } else if ( config->GetSmoothOnSurface() ) {

    for (unsigned long iVertex =0; iVertex<geometry->nVertex[val_marker]; iVertex++)  {

      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      for (iDim = 0; iDim < nDim; iDim++) {
        solver->GetNodes()->SetSensitivity(iPoint, iDim, LinSysSol[iVertex]);
      }
    }

  } else {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint*nDim + iDim;
        solver->GetNodes()->SetSensitivity(iPoint, iDim, LinSysSol[total_index]);
      }
    }

  }
}


void CGradientSmoothingSolver::MultiplyByVolumeDeformationStiffness(CGeometry *geometry, CSolver *solver, CVolumetricMovement *grid_movement, CConfig *config, bool Transpose) {

  LinSysRes.SetValZero();
  LinSysSol.SetValZero();

  // extract the sensitivities
  SetBoundaryDerivativesForMultiplication(geometry, solver, config, Transpose);

  for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
    for (auto iDim = 0; iDim < nDim ; iDim++) {
      LinSysRes.SetBlock(iPoint, iDim, 1.0);
    }
  }

  // get the matrix vector product class (implicite version of get matrix)
  CMatrixVectorProduct<su2double>* mat_vec = grid_movement->GetStiffnessMatrixVectorProduct(geometry, config, Transpose);

  // perform the matrix vector product
  (*mat_vec)(LinSysRes, LinSysSol);

  CSysMatrix<su2double>& StiffMatrix = grid_movement->GetStiffnessMatrix(geometry, config);
  ofstream matrix ("stiffness.dat");
  StiffMatrix.printMat(matrix);
  matrix.close();

  // set the calculated sensitivities
  Set_Sensitivities(geometry, solver, config);

  delete mat_vec;
}


void CGradientSmoothingSolver::SetBoundaryDerivativesForMultiplication(CGeometry *geometry, CSolver *solver, CConfig *config, bool Transpose) {
  unsigned short iDim, iMarker;
  unsigned long iPoint, total_index, iVertex;

  if ( Transpose ) {

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_SobolevBC(iMarker) == YES)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < nDim; iDim++) {
            total_index = iPoint*nDim + iDim;
            LinSysRes[total_index] = solver->GetNodes()->GetSensitivity(iPoint, iDim);
            LinSysSol[total_index] = solver->GetNodes()->GetSensitivity(iPoint, iDim);
          }
        }
      }
    }
    if (LinSysRes.norm() == 0.0) cout << "Warning: Derivatives are zero!" << endl;
  } else {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint*nDim + iDim;
        LinSysRes[total_index] = solver->GetNodes()->GetSensitivity(iPoint, iDim);
        LinSysSol[total_index] = solver->GetNodes()->GetSensitivity(iPoint, iDim);
      }
    }
  }
}


std::vector<std::vector<su2double>> CGradientSmoothingSolver::GetElementCoordinates(CGeometry *geometry, std::vector<unsigned long>& indexNode, int EL_KIND) {

  std::vector<std::vector<su2double>> Coord;

  switch (EL_KIND) {

  case EL_LINE:

    for(auto iNode=0; iNode<2; iNode++) {
      Coord.push_back(std::vector<su2double>());
      for(auto iDim=0; iDim<2; iDim++) {
        Coord[iNode].push_back( Get_ValCoord(geometry, indexNode[iNode], iDim) );
      }
    }
    break;

  case EL_TRIA:

    for(auto iNode=0; iNode<3; iNode++) {
      Coord.push_back(std::vector<su2double>());
      for(auto iDim=0; iDim<3; iDim++) {
        Coord[iNode].push_back( Get_ValCoord(geometry, indexNode[iNode], iDim) );
      }
    }
    break;

  case EL_TRIA2:

    for(auto iNode=0; iNode<3; iNode++) {
      Coord.push_back(std::vector<su2double>());
      for(auto iDim=0; iDim<3; iDim++) {
        Coord[iNode].push_back( Get_ValCoord(geometry, indexNode[iNode], iDim) );
      }
    }
    break;

  case EL_QUAD:

    for(auto iNode=0; iNode<4; iNode++) {
      Coord.push_back(std::vector<su2double>());
      for(auto iDim=0; iDim<3; iDim++) {
        Coord[iNode].push_back( Get_ValCoord(geometry, indexNode[iNode], iDim) );
      }
    }
    break;

  default:
    std::cout << "Type of element is not supported. " <<std::endl;

  }

  return Coord;

}


void CGradientSmoothingSolver::ApplyGradientSmoothingOnSurface(CGeometry *geometry, CSolver *solver, CNumerics **numerics, CConfig *config, unsigned long val_marker) {

  unsigned long iPoint;
  unsigned short iDim;

  /*--- Initialize vector and sparse matrix ---*/
  LinSysSol.Initialize(geometry->nVertex[val_marker], geometry->nVertex[val_marker], 1, 0.0);
  LinSysRes.Initialize(geometry->nVertex[val_marker], geometry->nVertex[val_marker], 1, 0.0);
  Jacobian.InitOwnConnectivity(geometry->nVertex[val_marker], 1, 1, val_marker, geometry, config);
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();

  auxVecInp.Initialize(geometry->nVertex[val_marker], geometry->nVertex[val_marker], 2, 0.0);
  auxVecRHS.Initialize(geometry->nVertex[val_marker], geometry->nVertex[val_marker], 1, 0.0);
  auxVecInp.SetValZero();
  auxVecRHS.SetValZero();

  /*--- here we can set a customized rhs for the solver if SobolevDebugMode is set ---*/
  for (unsigned long iVertex =0; iVertex<geometry->nVertex[val_marker]; iVertex++)  {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    su2double x = geometry->node[iPoint]->GetCoord(0);
    su2double y = geometry->node[iPoint]->GetCoord(1);
    for (iDim=0;iDim<nDim;iDim++) {
      auxVecInp.SetBlock(iVertex, iDim, sin(3.14159265359*x));
    }

    //iPoint, iDim, -2.0*3.14159265359*3.14159265359*sin(3.14159265359*x)*sin(3.14159265359*y) );
  }

  ofstream input ("input.txt");
  auxVecInp.printVec(input);
  input.close();

  Compute_Surface_StiffMatrix(geometry, numerics, config, val_marker);

  Compute_Surface_Residual(geometry, solver, config, val_marker);

  if ( config->GetDirichletSurfaceBound() ) {
    BC_Surface_Dirichlet(geometry, config, val_marker);
  }

  ofstream matrix ("matrix.dat");
  Jacobian.printMat(matrix);
  matrix.close();

  for (unsigned long iVertex =0; iVertex<geometry->nVertex[val_marker]; iVertex++)  {
    auxVecRHS.SetBlock(iVertex, 0, LinSysRes.GetBlock(iVertex,0));
  }

  ofstream rhs ("rhs.txt");
  auxVecRHS.printVec(rhs);
  rhs.close();

  Solve_Linear_System(geometry, config);

  Set_Sensitivities(geometry, solver, config, val_marker);

  ofstream output ("output.txt");
  LinSysSol.printVec(output);
  output.close();

}
