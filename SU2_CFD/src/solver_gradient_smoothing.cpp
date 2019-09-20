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
#include <algorithm>


CGradientSmoothingSolver::CGradientSmoothingSolver(CGeometry *geometry, CConfig *config) : CSolver(false,true) {

  unsigned short iDim, jDim;
  unsigned int iElem;

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

  /*--- Allocate element properties - only the index, to allow further integration with CFEASolver on a later stage ---*/
  element_properties = new CProperty*[nElement];
  for (iElem = 0; iElem < nElement; iElem++){
    element_properties[iElem] = new CProperty(iElem);
  }

  /*--- auxiliary subarrays and submatrices ---*/
  Residual = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Residual[iDim] = 0.0;
  Solution = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Solution[iDim] = 0.0;
  #ifndef CODI_FORWARD_TYPE
    passivedouble **matrixId;
  #else
    su2double **matrixId;
  #endif

  /*--- Matrices to impose boundary conditions ---*/

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

  /*--- Term ij of the Jacobian ---*/
  Jacobian_ij = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Jacobian_ij[iDim] = new su2double [nDim];
    for (jDim = 0; jDim < nDim; jDim++) {
      Jacobian_ij[iDim][jDim] = 0.0;
    }
  }

  /*--- linear system ---*/
  if ( config->GetSepDim() ) {
    LinSysSol.Initialize(nPoint, nPointDomain, 1, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, 1, 0.0);
    Jacobian.Initialize(nPoint, nPointDomain, 1, 1, false, geometry, config);
  } else {
    LinSysSol.Initialize(nPoint, nPointDomain, nDim, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nDim, 0.0);
    Jacobian.Initialize(nPoint, nPointDomain, nDim, nDim, false, geometry, config);
  }

}


CGradientSmoothingSolver::~CGradientSmoothingSolver(void) {

  unsigned short iDim;

  delete [] Residual;
  delete [] Solution;
  delete [] matrixId;

  if (element_container != NULL) {
    for (unsigned short iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (unsigned short jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        if (element_container[iVar][jVar] != NULL) delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] mZeros_Aux[iDim];
    delete [] mId_Aux[iDim];
    delete [] Jacobian_ij[iDim];
  }
  delete [] mZeros_Aux;
  delete [] mId_Aux;
  delete [] Jacobian_ij;

}


void CGradientSmoothingSolver::ApplyGradientSmoothing(CGeometry *geometry, CSolver *solver, CNumerics **numerics, CConfig *config) {

  dir = 0;

    CSysVector<su2double> auxVec;
    auxVec.Initialize(nPoint, nPointDomain, nDim, 0.0);

    auxVec.SetValZero();
    for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
        for(auto iDim=0; iDim<nDim; iDim++) {
        auxVec.SetBlock(iPoint, iDim, solver->node[iPoint]->GetSensitivity(iDim));
        }
    }

    ofstream input ("input.txt");
    auxVec.printVec(input);
    input.close();

  /*--- Initialize vector and sparse matrix ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();

  Compute_StiffMatrix(geometry, numerics, config);

  ofstream matrix ("Omatrix.dat");
  Jacobian.printMat(matrix);
  matrix.close();

  if ( config->GetSepDim() ) {

    for (dir = 0; dir < nDim ; dir++) {

      Compute_Residual(geometry, solver, config);
      Impose_BC(geometry, numerics, config);

      ofstream matrix ("matrix.dat");
      Jacobian.printMat(matrix);
      matrix.close();
      ofstream rhs ("rhs.txt");
      LinSysRes.printVec(rhs);
      rhs.close();

      Solve_Linear_System(geometry, config);
      Set_Sensitivities(geometry, solver, config);

      ofstream result ("result.txt");
      LinSysSol.printVec(result);
      result.close();

      LinSysSol.SetValZero();
      LinSysRes.SetValZero();
    }

  } else {

    Compute_Residual(geometry, solver, config);

    Impose_BC(geometry, numerics, config);

    ofstream rhs ("rhs.txt");
    LinSysRes.printVec(rhs);
    rhs.close();

    ofstream matrix ("matrix.dat");
    Jacobian.printMat(matrix);
    matrix.close();

    Solve_Linear_System(geometry, config);

    ofstream result ("result.txt");
    LinSysSol.printVec(result);
    result.close();

    Set_Sensitivities(geometry, solver, config);

  }

  auxVec.SetValZero();
  for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++)  {
      for(auto iDim=0; iDim<nDim; iDim++) {
      auxVec.SetBlock(iPoint, iDim, solver->node[iPoint]->GetSensitivity(iDim));
      }
  }

  ofstream output ("output.txt");
  auxVec.printVec(output);
  output.close();

  /*
  for (unsigned long iPoint =0; iPoint<geometry->GetnPoint(); iPoint++) {

    std::cout << "SU2 node "<< iPoint << " is mesh node " << geometry->node[iPoint]->GetGlobalIndex() << std::endl;

  }
  */
}


void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config){

  unsigned long iElem;
  unsigned short iNode, iDim, jDim, nNodes = 0;
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

        /*
         * std::cout << "Element: " << iElem
                  << " Jacobian Block: " << iNode << ", " << jNode << ", "
                  << indexNode[iNode] << ", " << indexNode[jNode] << ", "
                  << geometry->node[indexNode[iNode]]->GetGlobalIndex()<< ", " << geometry->node[indexNode[jNode]]->GetGlobalIndex() << std::endl;
        */

        if ( config->GetSepDim() ) {

          Jacobian_ij[0][0] = DHiDHj[dir][dir] + HiHj;
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij );

        } else {

          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0; jDim < nDim; jDim++) {
              Jacobian_ij[iDim][jDim] = DHiDHj[iDim][jDim];
              }
            Jacobian_ij[iDim][iDim] += HiHj;
          }
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

        }
      }

    }

  }

/*

  unsigned long iPoint;

  matrixId    = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++){
    matrixId[iDim]    = new su2double[nDim];
  }
  for(iDim = 0; iDim < nDim; iDim++){
    for (jDim = 0; jDim < nDim; jDim++){
      matrixId[iDim][jDim]    = 0.0;
    }
    matrixId[iDim][iDim] = 2.0;
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Jacobian.SetBlock(iPoint,iPoint,matrixId);
  }

*/
}


void CGradientSmoothingSolver::Compute_Residual(CGeometry *geometry, CSolver *solver, CConfig *config){

  unsigned long iElem;
  unsigned short iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
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

          Residual[dir] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * solver->node[indexNode[iNode]]->GetSensitivity(dir);
          LinSysRes.AddBlock(indexNode[iNode], &Residual[dir]);

        } else {

          for (iDim = 0; iDim < nDim; iDim++) {
            Residual[iDim] += Weight * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode,iGauss) * solver->node[indexNode[iNode]]->GetSensitivity(iDim);
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


void CGradientSmoothingSolver::Impose_BC(CGeometry *geometry, CNumerics **numerics, CConfig *config) {

  unsigned short iMarker;

  /*--- Get the boundary markers and iterate over them ---------------------------------*/
  /* Notice that for no marker we automatically impose Zero Neumann boundary conditions */

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_SobolevBC(iMarker) == YES) {
      BC_Dirichlet(geometry, NULL, numerics, config, iMarker);
      // BC_Neumann(geometry, NULL, numerics, config, iMarker);
    }

  }

}


void CGradientSmoothingSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, unsigned short val_marker) {


  unsigned long iPoint, iVertex;
  unsigned long iVar, jVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    // std::cout << "Dirichlet boundary in node: " << iPoint << " with marker: " << val_marker << " bound " << iVertex << std::endl;

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

// For now: left empty since there is no calculation necessary for zero Neumann boundaries.
void CGradientSmoothingSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, unsigned short val_marker) {

}


void CGradientSmoothingSolver::Solve_Linear_System(CGeometry *geometry, CConfig *config){

  unsigned long IterLinSol = 0;

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  SetIterLinSolver(IterLinSol);

  // valResidual = System.GetResidual();

}


void CGradientSmoothingSolver::Set_Sensitivities(CGeometry *geometry, CSolver *solver, CConfig *config){

  unsigned long iPoint, total_index;
  unsigned short iDim;

  if ( config->GetSepDim() ) {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      solver->node[iPoint]->SetSensitivity(dir,LinSysSol[iPoint]);

                std::cout << "Output " << iPoint << ", " << dir << ", " << LinSysSol[iPoint] << std::endl;

    }

  } else {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint*nDim + iDim;
        solver->node[iPoint]->SetSensitivity(iDim,LinSysSol[total_index]);
      }
    }
  }

}

