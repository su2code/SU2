/*!
 * \file CGradientSmoothing.cpp
 * \brief Main solver routines for the gradient smoothing problem.
 * \author T. Dick
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include <algorithm>

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CGradientSmoothingSolver.hpp"
#include "../../include/variables/CSobolevSmoothingVariable.hpp"

/*!
 * \namespace CGradientSmoothingSolverDetails
 * \brief Namespace to hide helper functions used by the CGradientSmoothingSolver class.
 */
namespace CGradientSmoothingSolverDetails {

template<typename su2matvecscalar>
void WriteVectorToGeometry(CGeometry* geometry, const CSysVector<su2matvecscalar>& vector) {
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++) {
    for (auto iDim = 0u; iDim < geometry->GetnDim(); iDim++) {
      geometry->SetSensitivity(iPoint,iDim, vector(iPoint, iDim));
    }
  }
}

template<typename su2matvecscalar>
void ReadVectorToGeometry(const CGeometry* geometry, CSysVector<su2matvecscalar>& vector) {
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++) {
    for (auto iDim = 0u; iDim < geometry->GetnDim(); iDim++) {
      vector(iPoint, iDim) = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint,iDim));
    }
  }
}
} // namespace CGradientSmoothingSolverDetails

CGradientSmoothingSolver::CGradientSmoothingSolver(CGeometry *geometry, CConfig *config) : CFEASolverBase(LINEAR_SOLVER_MODE::GRADIENT_MODE) {
  unsigned int marker_count = 0;
  unsigned long iPoint;

  /*--- This solver does not perform recording operations --*/
  adjoint = false;

  /*--- General geometric settings ---*/
  nDim         = geometry->GetnDim();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nElement     = geometry->GetnElem();

  /*--- Initialize the element container and assign the kind of each element ---*/

  element_container = new CElement** [MAX_TERMS]();
  for (unsigned int iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS]();

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

  /*--- For operations on surfaces we initalize the structures for nDim-1 ---*/
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

  /*--- initializations for linear equation systems ---*/
  if ( !config->GetSmoothOnSurface() ) {
    nVar = config->GetSmoothSepDim() ? 1 : nDim;
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config, false, true);
  } else {
    if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::PARAM_LEVEL_COMPLETE) {
      Jacobian.Initialize(nPoint, nPointDomain, nDim, nDim, false, geometry, config, false , true);
    } else {
      LinSysSol.Initialize(nPoint, nPointDomain, 1, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, 1, 0.0);
      Jacobian.Initialize(nPoint, nPointDomain, 1, 1, false, geometry, config, false, true);
    }
    visited.resize(geometry->GetnPoint(), false);
  }

  /*--- vectors needed for projection when working on the complete system ---*/
  if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::PARAM_LEVEL_COMPLETE ||
      config->GetSobMode() == ENUM_SOBOLEV_MODUS::ONLY_GRAD || config->GetSobMode() == ENUM_SOBOLEV_MODUS::MESH_LEVEL) {
    activeCoord.Initialize(nPoint, nPointDomain, nDim, 0.0);
    helperVecIn.Initialize(nPoint, nPointDomain, nDim, 0.0);
    helperVecOut.Initialize(nPoint, nPointDomain, nDim, 0.0);
    helperVecAux.Initialize(nPoint, nPointDomain, nDim, 0.0);
  }

  /*--- Initialize the CVariable structure holding solution data ---*/
  nodes = new CSobolevSmoothingVariable(nPoint, nDim, config);
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
  nodes->AllocateBoundaryVariables();

  /*--- vector for the parameter gradient ---*/
  deltaP.resize(config->GetnDV_Total(), 0.0);

  /*--- set the solver name for output purposes ---*/
  SolverName = "SOBOLEV";
}

CGradientSmoothingSolver::~CGradientSmoothingSolver() {

  delete nodes;

}

void CGradientSmoothingSolver::ApplyGradientSmoothingVolume(CGeometry* geometry, CNumerics* numerics,
                                                            const CConfig* config) {
  /*--- current dimension if we run consecutive on each dimension ---*/
  unsigned int iDim = 0;

  /*--- Set vector and sparse matrix to 0 ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();

  /*--- Compute the stiffness matrix for the smoothing operator. ---*/
  Compute_StiffMatrix(geometry, numerics, config);

  /*--- Impose boundary conditions to the RHS and solve the system. ---*/
  if (config->GetSmoothSepDim()) {
    for (iDim = 0; iDim < nDim ; iDim++) {

      SetCurrentDim(iDim);

      Compute_Residual(geometry, config);

      Impose_BC(geometry, config);

      Solve_Linear_System(geometry, config);

      WriteSensitivity(geometry, config);

      LinSysSol.SetValZero();
      LinSysRes.SetValZero();
    }

  } else {
    Compute_Residual(geometry, config);

    Impose_BC(geometry, config);

    Solve_Linear_System(geometry, config);

    WriteSensitivity(geometry, config);
  }
}

void CGradientSmoothingSolver::ApplyGradientSmoothingSurface(CGeometry* geometry, CNumerics* numerics,
                                                             const CConfig* config) {
  /*--- Set vector and sparse matrix to 0 ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();
  std::fill(visited.begin(), visited.end(), false);

  /*--- Loop over all DV markers to compute the stiffness matrix for the smoothing operator. ---*/
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      /*--- Compute the stiffness matrix for the smoothing operator. ---*/
      Compute_Surface_StiffMatrix(geometry, numerics, config, iMarker);

      Compute_Surface_Residual(geometry, config, iMarker);

      if (config->GetDirichletSurfaceBound()) {
        BC_Surface_Dirichlet(geometry, config, iMarker);
      }
    }
  }

  /*--- Set the matrix to identity if the current mpi rank holds no part of the DV marker. ---*/
  Complete_Surface_StiffMatrix(geometry);

  /*--- Solve the system and write the result back. ---*/
  Solve_Linear_System(geometry, config);

  WriteSensitivity(geometry, config);
}

void CGradientSmoothingSolver::ApplyGradientSmoothingDV(CGeometry* geometry, CNumerics* numerics,
                                                        CSurfaceMovement* surface_movement,
                                                        CVolumetricMovement* grid_movement, CConfig* config,
                                                        su2double** Gradient) {
  unsigned nDVtotal=config->GetnDV_Total();
  unsigned column, row;
  unsigned long iPoint;
  unsigned short iDim;
  vector<su2double> seedvector(nDVtotal, 0.0);
  vector<su2double> x(nDVtotal, 0.0);
  hessian.Initialize(nDVtotal);

  /*--- Reset the Jacobian to 0 ---*/
  Jacobian.SetValZero();

  /*--- Record the parameterization on the AD tape. ---*/
  if (rank == MASTER_NODE)  cout << " calculate the original gradient" << endl;
  RecordParameterizationJacobian(geometry, surface_movement, activeCoord, config);

  /*--- Calculate the original gradient. ---*/
  CalculateOriginalGradient(geometry, grid_movement, config, Gradient);

  /*--- Compute the full Sobolev Hessian approximation column by column. ---*/
  if (rank == MASTER_NODE)  cout << " computing the system matrix line by line" << endl;

  auto mat_vec = GetStiffnessMatrixVectorProduct<su2matvecscalar>(geometry, numerics, config);

  for (column=0; column<nDVtotal; column++) {

    if (rank == MASTER_NODE)  cout << "    working in column " << column << endl;

    /*--- Create a seeding vector and reset auxilliary vectors for the surface case. ---*/
    std::fill(seedvector.begin(), seedvector.end(), 0.0);
    seedvector[column] = 1.0;

    helperVecIn.SetValZero();
    helperVecOut.SetValZero();
    helperVecAux.SetValZero();

    /*--- Forward mode evaluation, i.e., matrix vector produt with the parameterization Jacobian.---*/
    ProjectDVtoMesh(geometry, seedvector, helperVecIn, activeCoord, config);

    /*--- Matrix vector product with the Laplace-Beltrami stiffness matrix. ---*/
    if (config->GetSmoothOnSurface()) {

      CSysMatrixComms::Initiate(helperVecIn, geometry, config, SOLUTION_MATRIX);
      CSysMatrixComms::Complete(helperVecIn, geometry, config, SOLUTION_MATRIX);

      mat_vec(helperVecIn, helperVecAux);

      /*--- Delete entries from halo cells, otherwise we get errors for parallel computations ---*/
      for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          helperVecOut(iPoint,iDim) = helperVecAux(iPoint,iDim);
        }
      }

    } else {

      /*--- Forward evaluation of the mesh deformation ---*/
      CGradientSmoothingSolverDetails::WriteVectorToGeometry<su2matvecscalar>(geometry, helperVecIn);
      grid_movement->SetVolume_Deformation(geometry, config, false, true, true);
      CGradientSmoothingSolverDetails::ReadVectorToGeometry<su2matvecscalar>(geometry, helperVecIn);

      CSysMatrixComms::Initiate(helperVecIn, geometry, config, SOLUTION_MATRIX);
      CSysMatrixComms::Complete(helperVecIn, geometry, config, SOLUTION_MATRIX);

      mat_vec(helperVecIn, helperVecAux);

      /*--- Only consider values inside the domain for the output of the matrix vector product ---*/
      for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          helperVecOut(iPoint,iDim) = helperVecAux(iPoint,iDim);
        }
      }

      /*--- Forward evaluation of the mesh deformation ---*/
      CGradientSmoothingSolverDetails::WriteVectorToGeometry<su2matvecscalar>(geometry, helperVecOut);
      grid_movement->SetVolume_Deformation(geometry, config, false, true, false);
      CGradientSmoothingSolverDetails::ReadVectorToGeometry<su2matvecscalar>(geometry, helperVecOut);

    }

    /*--- Reverse mode evaluation, i.e., matrix vector produt with the transposed parameterization Jacobian.---*/
    ProjectMeshToDV(geometry, helperVecOut, seedvector, activeCoord, config);

    /*--- Extract the projected direction ---*/
    for (row=0; row<nDVtotal; row++) {
      hessian(row,column) = SU2_TYPE::GetValue(seedvector[row]);
    }
  }

  /*--- Output the complete system matrix. ---*/
  if (rank == MASTER_NODE) {
    /*--- For multizone append zone number to filename. ---*/
    string hess_filename = config->GetObjFunc_Hess_FileName();
    if (config->GetMultizone_Problem()) {
      hess_filename = config->GetMultizone_FileName(hess_filename, config->GetiZone(), ".dat");
    }
    ofstream SysMatrix(hess_filename);
    SysMatrix.precision(config->GetOutput_Precision());

    for (row=0; row<nDVtotal; row++) {
      for (column=0; column<nDVtotal; column++) {
        SysMatrix << hessian(row,column);
        if (column!=nDVtotal-1) SysMatrix << ", ";
      }
      if (row!=nDVtotal-1) SysMatrix << '\n';
    }
    SysMatrix.close();
  }

  /*--- Calculate and output the treated gradient. ---*/
  hessian.Invert();
  hessian.MatVecMult(deltaP.begin(), x.begin());
  deltaP = x;

  OutputDVGradient(config);
}

void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry* geometry, CNumerics* numerics, const CConfig* config) {
  unsigned long iElem, iNode;
  unsigned int iDim, nNodes = 0, NelNodes, jNode;
  std::array<unsigned long, MAXNNODE_3D> indexNode;
  su2double val_Coord, HiHj = 0.0;
  int EL_KIND = 0;
  su2activematrix DHiDHj, Jacobian_block;

  Jacobian_block.resize(nDim, nDim) = su2double(0.0);

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

    numerics->Compute_Tangent_Matrix(element_container[GRAD_TERM][EL_KIND], config);

    NelNodes = element_container[GRAD_TERM][EL_KIND]->GetnNodes();

    /*--- for all nodes add the contribution to the system Jacobian ---*/

    for (iNode = 0; iNode < NelNodes; iNode++) {

      for (jNode = 0; jNode < NelNodes; jNode++) {

        DHiDHj = element_container[GRAD_TERM][EL_KIND]->Get_DHiDHj(iNode, jNode);
        HiHj = element_container[GRAD_TERM][EL_KIND]->Get_HiHj(iNode, jNode);

        if (config->GetSmoothSepDim()) {
          Jacobian_block[0][0] = DHiDHj[GetCurrentDim()][GetCurrentDim()] + HiHj;
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

void CGradientSmoothingSolver::Compute_Surface_StiffMatrix(CGeometry* geometry, CNumerics* numerics,
                                                           const CConfig* config, unsigned long val_marker,
                                                           unsigned int nSurfDim) {
  unsigned long iElem, iPoint, iVertex, iDim, iSurfDim;
  unsigned int iNode, jNode, nNodes = 0, NelNodes;
  std::array<unsigned long, MAXNNODE_2D> indexNode;
  std::array<unsigned long, MAXNNODE_2D> indexVertex;
  int EL_KIND = 0;
  su2double val_Coord, HiHj = 0.0;
  su2activematrix DHiDHj, Jacobian_block;

  Jacobian_block.resize(nDim, nDim) = su2double(0.0);

  /*--- Loops over all the elements ---*/
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    /*--- Identify the kind of boundary element ---*/
    GetElemKindAndNumNodes(geometry->bound[val_marker][iElem]->GetVTK_Type(), EL_KIND, nNodes);

    /*--- Retrieve the boundary reference and current coordinates ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element_container[GRAD_TERM][EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
      }
    }

    /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
      for (iNode = 0; iNode < nNodes; iNode++) {
        if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
      }
    }

    /*--- compute the contributions of the single elements inside the numerics container ---*/
    numerics->Compute_Tangent_Matrix(element_container[GRAD_TERM][EL_KIND], config);

    NelNodes = element_container[GRAD_TERM][EL_KIND]->GetnNodes();

    /*--- for all nodes add the contribution to the system Jacobian ---*/

    for (iNode = 0; iNode < NelNodes; iNode++) {
      for (jNode = 0; jNode < NelNodes; jNode++) {

        DHiDHj = element_container[GRAD_TERM][EL_KIND]->Get_DHiDHj(iNode, jNode);
        HiHj = element_container[GRAD_TERM][EL_KIND]->Get_HiHj(iNode, jNode);

        /*--- Get the calculated stiffness blocks for the current element and add them to the block.  ---*/
        for (iSurfDim=0; iSurfDim<nSurfDim; iSurfDim++) {
          Jacobian_block[iSurfDim][iSurfDim] = DHiDHj[0][0] + HiHj;
        }
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_block);

        /*--- Store if we set values for a given node already ---*/
        visited[indexNode[iNode]]=true;
        visited[indexNode[jNode]]=true;

      }
    }
  }
}

void CGradientSmoothingSolver::Compute_Residual(CGeometry* geometry, const CConfig* config) {
  unsigned long iElem;
  unsigned int iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  std::array<unsigned long, MAXNNODE_3D> indexNode;
  std::array<su2double, MAXNDIM> Residual;
  su2double Weight, Jac_X, val_Coord;

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
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
        if (config->GetSmoothSepDim()) {
          Residual[GetCurrentDim()] +=
              Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode, iGauss) * nodes->GetSensitivity(indexNode[iNode], GetCurrentDim());
          LinSysRes.AddBlock(indexNode[iNode], &Residual[GetCurrentDim()]);
        } else {
          for (iDim = 0; iDim < nDim; iDim++) {
            Residual[iDim] += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode, iGauss) * nodes->GetSensitivity(indexNode[iNode], iDim);
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

void CGradientSmoothingSolver::Compute_Surface_Residual(CGeometry* geometry, const CConfig* config,
                                                        unsigned long val_marker) {
  unsigned long iElem, iPoint, iVertex;
  unsigned int iDim, iNode, nNodes = 0;
  int EL_KIND = 0;
  std::array<unsigned long, MAXNNODE_2D> indexNode;
  std::array<unsigned long, MAXNNODE_2D> indexVertex;
  su2double Weight, Jac_X, norm, val_Coord, normalSens = 0.0, Residual=0.0;
  su2double normal[MAXNDIM];

  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    /*--- Identify the kind of boundary element ---*/
    GetElemKindAndNumNodes(geometry->bound[val_marker][iElem]->GetVTK_Type(), EL_KIND, nNodes);

    /*--- Retrieve the boundary reference and current coordinates ---*/
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim);
        element_container[GRAD_TERM][EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
      }
    }

    /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
      for (iNode = 0; iNode < nNodes; iNode++) {
        if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
      }
    }

    element_container[GRAD_TERM][EL_KIND]
        ->ClearElement(); /*--- Restarts the element: avoids adding over previous results in other elements --*/
    element_container[GRAD_TERM][EL_KIND]->ComputeGrad_SurfaceEmbedded();
    unsigned int nGauss = element_container[GRAD_TERM][EL_KIND]->GetnGaussPoints();

    for (unsigned int iGauss = 0; iGauss < nGauss; iGauss++) {
      Weight = element_container[GRAD_TERM][EL_KIND]->GetWeight(iGauss);
      Jac_X = element_container[GRAD_TERM][EL_KIND]->GetJ_X(iGauss);

      for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
        geometry->vertex[val_marker][indexVertex[iNode]]->GetNormal(normal);
        norm = GeometryToolbox::Norm(nDim, normal);
        for (iDim = 0; iDim < nDim; iDim++) {
          normal[iDim] = normal[iDim] / norm;
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          normalSens += normal[iDim] * nodes->GetSensitivity(indexNode[iNode], iDim);
        }

        Residual += Weight * Jac_X * element_container[GRAD_TERM][EL_KIND]->GetNi(iNode, iGauss) * normalSens;
        LinSysRes.AddBlock(indexNode[iNode], &Residual);

        Residual = 0;
        normalSens = 0;
      }
    }
  }
}

void CGradientSmoothingSolver::Impose_BC(const CGeometry* geometry, const CConfig* config) {
  unsigned int iMarker;

  /*--- Get the boundary markers and iterate over them
   * This means for no specified marker we automatically impose Zero Neumann boundary conditions. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_SobolevBC(iMarker) == YES) {
      BC_Dirichlet(geometry, config, iMarker);
    }
  }
}

void CGradientSmoothingSolver::BC_Dirichlet(const CGeometry* geometry, const CConfig* config, unsigned int val_marker) {
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

void CGradientSmoothingSolver::BC_Surface_Dirichlet(const CGeometry* geometry, const CConfig* config,
                                                    unsigned int val_marker) {
  unsigned long iPoint, iVertex;
  const su2double zeros[MAXNDIM] = {0.0};

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (nodes->GetIsBoundaryPoint(iPoint)) {
      LinSysSol.SetBlock(iPoint, zeros);
      LinSysRes.SetBlock(iPoint, zeros);
      Jacobian.EnforceSolutionAtNode(iPoint, zeros, LinSysRes);
    }
  }
}

void CGradientSmoothingSolver::Solve_Linear_System(CGeometry* geometry, const CConfig* config) {
  /* For MPI prescribe vector entries across the ranks before solving the system.
   * Analog to FEA solver this is only done for the solution */
  CSysMatrixComms::Initiate(LinSysSol, geometry, config);
  CSysMatrixComms::Complete(LinSysSol, geometry, config);

  /*--- elimination shedule ---*/
  Set_VertexEliminationSchedule(geometry, config);

  SU2_OMP_PARALLEL
  {

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  SU2_OMP_MASTER
  {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  END_SU2_OMP_MASTER
  }
  END_SU2_OMP_PARALLEL
}

template <typename scalar_type>
CSysMatrixVectorProduct<scalar_type> CGradientSmoothingSolver::GetStiffnessMatrixVectorProduct(CGeometry* geometry,
                                                                                               CNumerics* numerics,
                                                                                               const CConfig* config) {

  /*--- Compute the sparse stiffness matrix ---*/
  if (config->GetSmoothOnSurface()) {
    for (unsigned int iMarker = 0; iMarker < config->GetnMarker_CfgFile(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Compute_Surface_StiffMatrix(geometry, numerics, config, iMarker, nDim);
      }
    }
    Complete_Surface_StiffMatrix(geometry);
  } else {
    Compute_StiffMatrix(geometry, numerics, config);
  }

  return CSysMatrixVectorProduct<scalar_type>(Jacobian, geometry, config);
}

void CGradientSmoothingSolver::CalculateOriginalGradient(CGeometry *geometry, CVolumetricMovement *grid_movement, CConfig *config, su2double** Gradient) {

  unsigned int iDV, iDV_Value, iDV_index;

  if (rank == MASTER_NODE) cout << endl << "Calculating the original DV gradient." << endl;

  WriteSensToGeometry(geometry);

  grid_movement->SetVolume_Deformation(geometry, config, false, true);

  CGradientSmoothingSolverDetails::ReadVectorToGeometry<su2matvecscalar>(geometry, helperVecOut);

  ProjectMeshToDV(geometry, helperVecOut, deltaP, activeCoord, config);

  OutputDVGradient(config, "orig_grad.dat");

  iDV_index = 0;
  for (iDV = 0; iDV  < config->GetnDV(); iDV++){
    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++){
      Gradient[iDV][iDV_Value] = deltaP[iDV_index];
      iDV_index++;
    }
  }
}

void CGradientSmoothingSolver::RecordTapeAndCalculateOriginalGradient(CGeometry *geometry, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement, CConfig *config, su2double **Gradient) {

  /*--- Record the parameterization on the AD tape. ---*/
  if (rank == MASTER_NODE)  cout << " calculate the original gradient" << endl;
  RecordParameterizationJacobian(geometry, surface_movement, activeCoord, config);

  /*--- Calculate the original gradient. ---*/
  CalculateOriginalGradient(geometry, grid_movement, config, Gradient);

}

void CGradientSmoothingSolver::OutputDVGradient(const CConfig* config, string out_file) {
  unsigned iDV;
  if (rank == MASTER_NODE) {
    /*--- For multizone append zone number to filename. ---*/
    if (config->GetMultizone_Problem()) {
      out_file = config->GetMultizone_FileName(out_file, config->GetiZone(), ".dat");
    }
    ofstream delta_p (out_file);
    delta_p.precision(config->GetOutput_Precision());

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

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nPoint  = geometry->GetnPoint();
  nDV     = config->GetnDV();

  /*--- Start recording of operations ---*/

  AD::Reset();

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero,
   * since we want to have the derivative at alpha = 0, i.e. for the current design. ---*/

  for (iDV = 0; iDV < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization of su2double with 0.0 resets the existing AD index. ---*/
      config->SetDV_Value(iDV, iDV_Value, 0.0);
      AD::RegisterInput(config->GetDV_Value(iDV, iDV_Value));
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

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  /*--- Seeding for the DV. --*/
  iDV_index = 0;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      SU2_TYPE::SetDerivative(config->GetDV_Value(iDV, iDV_Value), SU2_TYPE::GetValue(seeding[iDV_index]));
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

  /*--- Compute derivatives and extract the gradient ---*/
  AD::ComputeAdjoint();

  iDV_index = 0;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      my_Gradient = SU2_TYPE::GetDerivative(config->GetDV_Value(iDV, iDV_Value));

      SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

      output[iDV_index] = localGradient;
      iDV_index++;
    }
  }

  AD::ClearAdjoints();

}

void CGradientSmoothingSolver::WriteSensToGeometry(CGeometry* geometry) const {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->SetSensitivity(iPoint, iDim, nodes->GetSensitivity(iPoint, iDim));
    }
  }
}

void CGradientSmoothingSolver::ReadSensFromGeometry(const CGeometry* geometry) {
  unsigned long iPoint;
  unsigned int iDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      nodes->SetSensitivity(iPoint, iDim, geometry->GetSensitivity(iPoint, iDim));
    }
  }
}

void CGradientSmoothingSolver::WriteSensitivity(CGeometry* geometry, const CConfig* config) {

  unsigned int iDim;
  su2double normal[MAXNDIM];

  /*--- Split between surface and volume first. ---*/
  if ( config->GetSmoothOnSurface() ) {
    /*--- Write back values for all design boundaries in the surface case. ---*/
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES) {
        for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->GetNormal(normal);
          const auto norm = GeometryToolbox::Norm(nDim, normal);
          for (iDim = 0; iDim < nDim; iDim++) {
            normal[iDim] = normal[iDim] / norm;
            nodes->SetSensitivity(iPoint, iDim, normal[iDim] * LinSysSol[iPoint]);
          }
        }
      }
    }
  } else {
    if (config->GetSmoothSepDim()) {
      for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        nodes->SetSensitivity(iPoint, GetCurrentDim(), LinSysSol[iPoint]);
      }
    } else {
      for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          nodes->SetSensitivity(iPoint, iDim, LinSysSol(iPoint, iDim));
        }
      }
    }
  }
}

void CGradientSmoothingSolver::Set_VertexEliminationSchedule(CGeometry* geometry, const CConfig* config) {
  /*--- Store global point indices of essential BC markers. ---*/
  vector<unsigned long> myPoints;

  /*--- get global indexes of Dirichlet boundaries ---*/

  if (config->GetSmoothOnSurface()) {
    /*--- Surface case:
     * Fix design boundary border if Dirichlet condition is set for the DV marker
     * and the current rank holds part of it. ---*/
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES && config->GetDirichletSurfaceBound()) {
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          /*--- Get node index ---*/
          const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (nodes->GetIsBoundaryPoint(iPoint)) {
            myPoints.push_back(geometry->nodes->GetGlobalIndex(iPoint));
          }
        }
      }
    }

  } else {

    /*--- Volume case:
     * All boundaries where Dirichlet conditions are set. ---*/
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

  CommunicateExtraEliminationVertices(geometry, myPoints);

  /*--- eliminate the extra vertices ---*/
  for (auto iPoint : ExtraVerticesToEliminate) {
    Jacobian.EnforceSolutionAtNode(iPoint, LinSysSol.GetBlock(iPoint), LinSysRes);
  }
}

void CGradientSmoothingSolver::Complete_Surface_StiffMatrix(const CGeometry* geometry) {

  /*--- Assembling the stiffness matrix on the design surface means the Jacobian is the identity for nodes inside the domain. ---*/
  for (unsigned long iPoint = 0ul; iPoint < geometry->GetnPointDomain(); iPoint++){
    if (!visited[iPoint]) {
      Jacobian.AddVal2Diag(iPoint, 1.0);
    }
  }

}
