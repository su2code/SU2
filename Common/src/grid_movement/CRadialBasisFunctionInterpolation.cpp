/*!
 * \file CRadialBasisFunctionInterpolation.cpp
 * \brief Subroutines for moving mesh volume elements using Radial Basis Function interpolation.
 * \author F. van Steen
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/grid_movement/CRadialBasisFunctionInterpolation.hpp"
#include "../../include/interface_interpolation/CRadialBasisFunction.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"
#include "../../include/toolboxes/CSymmetricMatrix.hpp"

CRadialBasisFunctionInterpolation::CRadialBasisFunctionInterpolation(CGeometry* geometry, CConfig* config)
    : CVolumetricMovement(geometry) {}

CRadialBasisFunctionInterpolation::~CRadialBasisFunctionInterpolation() {
  for (auto*& ptr : BoundNodes) {
    delete ptr;
    ptr = nullptr;
  }
}

void CRadialBasisFunctionInterpolation::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo,
                                                              bool Derivative, bool ForwardProjectionDerivative) {
  /*--- Retrieve type of RBF and its support radius ---*/

  const auto kindRBF = config->GetKindRadialBasisFunction();
  const su2double radius = config->GetRadialBasisFunctionParameter();

  su2double MinVolume, MaxVolume;

  /*--- Retrieving number of deformation steps and screen output from config ---*/

  const auto Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();
  auto Screen_Output = config->GetDeform_Output();

  /*--- Disable the screen output if we're running SU2_CFD ---*/

  if (config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD && !Derivative) Screen_Output = false;
  if (config->GetSmoothGradient()) Screen_Output = true;

  /*--- Determining the boundary and internal nodes. Setting the control nodes. ---*/
  SetBoundNodes(geometry, config);

  vector<unsigned long> internalNodes;
  SetInternalNodes(geometry, config, internalNodes);

  SetCtrlNodes(config);

  /*--- Looping over the number of deformation iterations ---*/
  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    /*--- Compute min volume in the entire mesh. ---*/

    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
    if (rank == MASTER_NODE && Screen_Output)
      cout << "Min. volume: " << MinVolume << ", max. volume: " << MaxVolume << "." << endl;

    /*--- Solving the RBF system, resulting in the interpolation coefficients ---*/
    SolveRBFSystem(geometry, config, kindRBF, radius);

    /*--- Updating the coordinates of the grid ---*/
    UpdateGridCoord(geometry, config, kindRBF, radius, internalNodes);

    if (UpdateGeo) {
      UpdateDualGrid(geometry, config);
    }

    /*--- Check for failed deformation (negative volumes). ---*/

    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);

    /*--- Calculate amount of nonconvex elements ---*/

    ComputenNonconvexElements(geometry, Screen_Output);

    if (rank == MASTER_NODE && Screen_Output) {
      cout << "Non-linear iter.: " << iNonlinear_Iter + 1 << "/" << Nonlinear_Iter << ". ";
      if (nDim == 2)
        cout << "Min. area: " << MinVolume << "." << endl;
      else
        cout << "Min. volume: " << MinVolume << "." << endl;
    }
  }
}

void CRadialBasisFunctionInterpolation::SolveRBFSystem(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type,
                                                       const su2double radius) {
  /*--- In case of data reduction an iterative greedy algorithm is applied
          to perform the interpolation with a reduced set of control nodes.
          Otherwise with a full set of control nodes. ---*/
  const auto& rbfParam = config->GetRBFParam();

  if (rbfParam.DataReduction) {
    /*--- Error tolerance for the data reduction tolerance ---*/
    const su2double dataReductionTolerance = rbfParam.GreedyTolerance;

    /*--- Local maximum error node and corresponding maximum error  ---*/
    unsigned long maxErrorNodeLocal;
    su2double maxErrorLocal{0};

    /*--- Obtaining the initial maximum error nodes, which are found based on the maximum applied deformation. */
    if (ControlNodes->empty()) {
      GetInitMaxErrorNode(geometry, config, maxErrorNodeLocal, maxErrorLocal);
      SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    }

    /*--- Number of greedy iterations. ---*/
    unsigned short greedyIter = 0;

    /*--- While the maximum error is above the tolerance, data reduction algorithm is continued. ---*/
    while (MaxErrorGlobal > dataReductionTolerance || greedyIter == 0) {
      /*--- In case of a nonzero local error, control nodes are added ---*/
      if (maxErrorLocal > 0) {
        AddControlNode(maxErrorNodeLocal);
      }

      /*--- Obtaining the global number of control nodes. ---*/
      ComputeNCtrlNodesGlobal();

      /*--- Obtaining the interpolation coefficients. ---*/
      GetInterpCoeffs(geometry, config, type, radius);

      /*--- Determining the interpolation error, of the non-control boundary nodes. ---*/
      GetInterpError(geometry, config, type, radius, maxErrorNodeLocal, maxErrorLocal);
      SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

      if (rank == MASTER_NODE) {
        cout << "Greedy iteration: " << greedyIter << ". Max error: " << MaxErrorGlobal
             << ". Global nr. of ctrl nodes: " << nCtrlNodesGlobal << "\n"
             << endl;
      }
      greedyIter++;
    }
  } else {
    /*--- Obtaining the interpolation coefficients. ---*/
    GetInterpCoeffs(geometry, config, type, radius);
  }
}

void CRadialBasisFunctionInterpolation::GetInterpCoeffs(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type,
                                                        const su2double radius) {
  /*--- Obtaining the control nodes coordinates and distributing over all processes. ---*/
  SetCtrlNodeCoords(geometry);

  /*--- Obtaining the deformation of the control nodes. ---*/
  SetDeformation(geometry, config);

  /*--- Computation of the (inverse) interpolation matrix. ---*/
  su2passivematrix invInterpMat;
  ComputeInterpolationMatrix(geometry, type, radius, invInterpMat);

  /*--- Obtaining the interpolation coefficients. ---*/
  ComputeInterpCoeffs(invInterpMat);
}

void CRadialBasisFunctionInterpolation::SetBoundNodes(CGeometry* geometry, CConfig* config) {
  /*--- Storing of the local node, marker and vertex information of the boundary nodes ---*/

  /*--- Looping over the markers ---*/
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    /*--- Checking if not internal or send/receive marker ---*/
    if (!config->GetMarker_All_Deform_Mesh_Internal(iMarker) && !config->GetMarker_All_SendRecv(iMarker)) {
      /*--- Looping over the vertices of marker ---*/
      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        /*--- Node in consideration ---*/
        auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Check whether node is part of the subdomain and not shared with a receiving marker (for parallel
         * computation) ---*/
        if (geometry->nodes->GetDomain(iNode)) {
          BoundNodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
        }
      }
    }
  }

  /*--- Sorting of the boundary nodes based on their index ---*/
  const auto smaller_index = [](const CRadialBasisFunctionNode* a, const CRadialBasisFunctionNode* b) {
    return a->GetIndex() < b->GetIndex();
  };
  sort(BoundNodes.begin(), BoundNodes.end(), smaller_index);

  /*--- Obtaining unique set ---*/
  const auto equal_index = [](const CRadialBasisFunctionNode* a, const CRadialBasisFunctionNode* b) {
    return a->GetIndex() == b->GetIndex();
  };
  BoundNodes.resize(std::distance(BoundNodes.begin(), unique(BoundNodes.begin(), BoundNodes.end(), equal_index)));
}

void CRadialBasisFunctionInterpolation::SetCtrlNodes(CConfig* config) {
  /*--- Assigning the control nodes based on whether data reduction is applied or not. ---*/
  if (config->GetRBFParam().DataReduction) {
    /*--- Control nodes are an empty set ---*/
    ControlNodes = &ReducedControlNodes;
  } else {
    /*--- Control nodes are the boundary nodes ---*/
    ControlNodes = &BoundNodes;
  }

  /*--- Obtaining the total number of control nodes. ---*/
  ComputeNCtrlNodesGlobal();
}

void CRadialBasisFunctionInterpolation::ComputeInterpolationMatrix(CGeometry* geometry, const RADIAL_BASIS& type,
                                                                   const su2double radius,
                                                                   su2passivematrix& invInterpMat) {
  /*--- In case of parallel computation, the interpolation coefficients are computed on the master node ---*/

  if (rank == MASTER_NODE) {
    CSymmetricMatrix interpMat;

    /*--- Initialization of the interpolation matrix ---*/
    interpMat.Initialize(nCtrlNodesGlobal);

    /*--- Construction of the interpolation matrix.
      Since this matrix is symmetric only upper halve has to be considered ---*/

    /*--- Looping over the target nodes ---*/
    for (auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++) {
      /*--- Looping over the control nodes ---*/
      for (auto jNode = iNode; jNode < nCtrlNodesGlobal; jNode++) {
        /*--- Distance between nodes ---*/
        auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[iNode * nDim], CtrlCoords[jNode * nDim]);

        /*--- Evaluation of RBF ---*/
        interpMat(iNode, jNode) = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
      }
    }

    /*--- Obtaining lower halve using symmetry ---*/
    const bool kernelIsSPD = (type == RADIAL_BASIS::WENDLAND_C2) || (type == RADIAL_BASIS::GAUSSIAN) ||
                             (type == RADIAL_BASIS::INV_MULTI_QUADRIC);

    /*--- inverting the interpolation matrix ---*/
    interpMat.Invert(kernelIsSPD);
    invInterpMat = interpMat.StealData();
  }
}

void CRadialBasisFunctionInterpolation::SetDeformation(CGeometry* geometry, CConfig* config) {
  /* --- Initialization of the deformation vector ---*/
  CtrlNodeDeformation.resize(ControlNodes->size() * nDim, 0.0);

  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  /*--- Loop over the control nodes ---*/
  for (auto iNode = 0ul; iNode < ControlNodes->size(); iNode++) {
    const auto iMarker = (*ControlNodes)[iNode]->GetMarker();

    /*--- Setting nonzero displacement of the deformation markers, else setting zero displacement for static markers
     * ---*/
    if (IsDeformationMarker(config, iMarker)) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        CtrlNodeDeformation[iNode * nDim + iDim] = SU2_TYPE::GetValue(
            geometry->vertex[iMarker][(*ControlNodes)[iNode]->GetVertex()]->GetVarCoord()[iDim] * VarIncrement);
      }
    } else {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        CtrlNodeDeformation[iNode * nDim + iDim] = 0.0;
      }
    }
  }

/*--- In case of a parallel computation, the deformation of all control nodes is send to the master process ---*/
#ifdef HAVE_MPI

  /*--- Local number of control nodes ---*/
  unsigned long Local_nControlNodes = ControlNodes->size();

  /*--- Array containing the local number of control nodes ---*/
  std::vector<unsigned long> Local_nControlNodesArr(size);

  /*--- gathering local control node coordinate sizes on all processes. ---*/
  SU2_MPI::Allgather(&Local_nControlNodes, 1, MPI_UNSIGNED_LONG, Local_nControlNodesArr.data(), 1, MPI_UNSIGNED_LONG,
                     SU2_MPI::GetComm());

  /*--- Gathering all deformation vectors on the master node ---*/
  if (rank == MASTER_NODE) {
    /*--- resizing the global deformation vector ---*/
    CtrlNodeDeformation.resize(nCtrlNodesGlobal * nDim);

    /*--- Receiving the local deformation vector from other processes ---*/
    unsigned long start_idx = 0;
    for (auto iProc = 0; iProc < size; iProc++) {
      if (iProc != MASTER_NODE) {
        SU2_MPI::Recv(CtrlNodeDeformation.data() + start_idx, Local_nControlNodesArr[iProc] * nDim, MPI_DOUBLE, iProc,
                      0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      }
      start_idx += Local_nControlNodesArr[iProc] * nDim;
    }

  } else {
    /*--- Sending the local deformation vector to the master node ---*/
    SU2_MPI::Send(CtrlNodeDeformation.data(), Local_nControlNodes * nDim, MPI_DOUBLE, MASTER_NODE, 0,
                  SU2_MPI::GetComm());
  }
#endif
}

void CRadialBasisFunctionInterpolation::SetInternalNodes(CGeometry* geometry, CConfig* config,
                                                         vector<unsigned long>& internalNodes) {
  /*--- Looping over all nodes and check if part of domain and not on boundary ---*/
  for (auto iNode = 0ul; iNode < geometry->GetnPoint(); iNode++) {
    if (!geometry->nodes->GetBoundary(iNode)) {
      internalNodes.push_back(iNode);
    }
  }

  /*--- Adding nodes on markers considered as internal nodes ---*/
  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++) {
    /*--- Check if marker is considered as internal nodes ---*/
    if (config->GetMarker_All_Deform_Mesh_Internal(iMarker)) {
      /*--- Loop over marker vertices ---*/
      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        /*--- Local node index ---*/
        auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- if not among the boundary nodes ---*/
        if (find_if(BoundNodes.begin(), BoundNodes.end(),
                    [&](CRadialBasisFunctionNode* i) { return i->GetIndex() == iNode; }) == BoundNodes.end()) {
          internalNodes.push_back(iNode);
        }
      }
    }
  }

  /*--- In case of a parallel computation, the nodes on the send/receive markers are included as internal nodes
          if they are not already a boundary node with known deformation ---*/

#ifdef HAVE_MPI
  /*--- Looping over the markers ---*/
  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++) {
    /*--- If send or receive marker ---*/
    if (config->GetMarker_All_SendRecv(iMarker)) {
      /*--- Loop over marker vertices ---*/
      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        /*--- Local node index ---*/
        auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- if not among the boundary nodes ---*/
        if (find_if(BoundNodes.begin(), BoundNodes.end(),
                    [&](CRadialBasisFunctionNode* i) { return i->GetIndex() == iNode; }) == BoundNodes.end()) {
          internalNodes.push_back(iNode);
        }
      }
    }
  }

  /*--- sorting of the local indices ---*/
  sort(internalNodes.begin(), internalNodes.end());

  /*--- Obtaining unique set of internal nodes ---*/
  internalNodes.resize(std::distance(internalNodes.begin(), unique(internalNodes.begin(), internalNodes.end())));
#endif
}

void CRadialBasisFunctionInterpolation::ComputeInterpCoeffs(su2passivematrix& invInterpMat) {
  /*--- resizing the interpolation coefficient vector ---*/
  InterpCoeff.resize(nDim * nCtrlNodesGlobal);

  /*--- Coefficients are found on the master process.
          Resulting coefficient is found by summing the multiplications of inverse interpolation matrix entries with
     deformation ---*/
  if (rank == MASTER_NODE) {
    for (auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        InterpCoeff[iNode * nDim + iDim] = 0;
        for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++) {
          InterpCoeff[iNode * nDim + iDim] += invInterpMat(iNode, jNode) * CtrlNodeDeformation[jNode * nDim + iDim];
        }
      }
    }
  }

/*--- Broadcasting the interpolation coefficients ---*/
#ifdef HAVE_MPI
  SU2_MPI::Bcast(InterpCoeff.data(), InterpCoeff.size(), MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
#endif
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type,
                                                        const su2double radius,
                                                        const vector<unsigned long>& internalNodes) {
  if (rank == MASTER_NODE) {
    cout << "updating the grid coordinates" << endl;
  }

  /*--- Update of internal node coordinates ---*/
  UpdateInternalCoords(geometry, type, radius, internalNodes);

  /*--- Update of boundary node coordinates ---*/
  UpdateBoundCoords(geometry, config, type, radius);

  /*--- In case of data reduction, perform the correction for nonzero error nodes ---*/
  if (config->GetRBFParam().DataReduction && BoundNodes.size() > 0) {
    SetCorrection(geometry, config, type, internalNodes);
  }
}

void CRadialBasisFunctionInterpolation::UpdateInternalCoords(CGeometry* geometry, const RADIAL_BASIS& type,
                                                             const su2double radius,
                                                             const vector<unsigned long>& internalNodes) {
  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[3] = {0.0};

  /*--- Loop over the internal nodes ---*/
  for (auto iNode = 0ul; iNode < internalNodes.size(); iNode++) {
    /*--- Loop for contribution of each control node ---*/
    for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++) {
      /*--- Determine distance between considered internal and control node ---*/
      auto dist =
          GeometryToolbox::Distance(nDim, CtrlCoords[jNode * nDim], geometry->nodes->GetCoord(internalNodes[iNode]));

      /*--- Evaluate RBF based on distance ---*/
      auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

      /*--- Add contribution to total coordinate variation ---*/
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        var_coord[iDim] += rbf * InterpCoeff[jNode * nDim + iDim];
      }
    }

    /*--- Apply the coordinate variation and resetting the var_coord vector to zero ---*/
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, var_coord[iDim]);
      var_coord[iDim] = 0;
    }
  }
}

void CRadialBasisFunctionInterpolation::UpdateBoundCoords(CGeometry* geometry, CConfig* config,
                                                          const RADIAL_BASIS& type, const su2double radius) {
  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[3] = {0.0};

  /*--- In case of data reduction, the non-control boundary nodes are treated as if they where internal nodes ---*/
  if (config->GetRBFParam().DataReduction) {
    /*--- Looping over the non selected boundary nodes ---*/
    for (auto iNode = 0ul; iNode < BoundNodes.size(); iNode++) {
      /*--- Finding contribution of each control node ---*/
      for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++) {
        /*--- Distance of non-selected boundary node to control node ---*/
        auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[jNode * nDim],
                                              geometry->nodes->GetCoord(BoundNodes[iNode]->GetIndex()));

        /*--- Evaluation of the radial basis function based on the distance ---*/
        auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

        /*--- Computing and add the resulting coordinate variation ---*/
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          var_coord[iDim] += rbf * InterpCoeff[jNode * nDim + iDim];
        }
      }

      /*--- Applying the coordinate variation and resetting the var_coord vector*/
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        geometry->nodes->AddCoord(BoundNodes[iNode]->GetIndex(), iDim, var_coord[iDim]);
        var_coord[iDim] = 0;
      }
    }
  }

  /*--- Applying the surface deformation, which are stored in the deformation vector ---*/
  for (auto jNode = 0ul; jNode < ControlNodes->size(); jNode++) {
    if (config->GetMarker_All_Moving((*ControlNodes)[jNode]->GetMarker())) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        geometry->nodes->AddCoord((*ControlNodes)[jNode]->GetIndex(), iDim, CtrlNodeDeformation[jNode * nDim + iDim]);
      }
    }
  }
}

void CRadialBasisFunctionInterpolation::GetInitMaxErrorNode(CGeometry* geometry, CConfig* config,
                                                            unsigned long& maxErrorNodeLocal,
                                                            su2double& maxErrorLocal) {
  /*--- Set max error to zero ---*/
  maxErrorLocal = 0.0;

  /*--- Loop over the nodes ---*/
  for (auto iNode = 0ul; iNode < BoundNodes.size(); iNode++) {
    /*--- Compute to squared norm of the deformation ---*/
    su2double normSquaredDeformation = GeometryToolbox::SquaredNorm(
        nDim, geometry->vertex[BoundNodes[iNode]->GetMarker()][BoundNodes[iNode]->GetVertex()]->GetVarCoord());

    /*--- In case squared norm deformation is larger than the error, update the error ---*/
    if (normSquaredDeformation > maxErrorLocal) {
      maxErrorLocal = normSquaredDeformation;
      maxErrorNodeLocal = iNode;
    }
  }

  /*--- Account for the possibility of applying the deformation in multiple steps ---*/
  maxErrorLocal = sqrt(maxErrorLocal) / ((su2double)config->GetGridDef_Nonlinear_Iter());
}

void CRadialBasisFunctionInterpolation::SetCtrlNodeCoords(CGeometry* geometry) {
  /*--- The coordinates of all control nodes are made available on all processes ---*/

  /*--- resizing the matrix containing the global control node coordinates ---*/
  std::cout << nCtrlNodesGlobal << std::endl;
  CtrlCoords.resize(nCtrlNodesGlobal * nDim);

  /*--- Array containing the local control node coordinates ---*/
  std::vector<su2double> localCoords(nDim * ControlNodes->size());

  /*--- Storing local control node coordinates ---*/
  for (auto iNode = 0ul; iNode < ControlNodes->size(); iNode++) {
    auto coord = geometry->nodes->GetCoord((*ControlNodes)[iNode]->GetIndex());
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      localCoords[iNode * nDim + iDim] = coord[iDim];
    }
  }

  /*--- Gathering local control node coordinate sizes on all processes. ---*/
  std::vector<int> LocalCoordsSizes(size);
  int localCoordsSize = nDim * ControlNodes->size();
  SU2_MPI::Allgather(&localCoordsSize, 1, MPI_INT, LocalCoordsSizes.data(), 1, MPI_INT, SU2_MPI::GetComm());

  /*--- Array containing the starting indices for the allgatherv operation */
  std::vector<int> disps(size, 0);

  for (auto iProc = 1; iProc < SU2_MPI::GetSize(); iProc++) {
    disps[iProc] = disps[iProc - 1] + LocalCoordsSizes[iProc - 1];
  }

  /*--- Distributing global control node coordinates among all processes ---*/
  SU2_MPI::Allgatherv(&localCoords, localCoordsSize, MPI_DOUBLE, CtrlCoords.data(), LocalCoordsSizes.data(),
                      disps.data(), MPI_DOUBLE, SU2_MPI::GetComm());
}

void CRadialBasisFunctionInterpolation::GetInterpError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type,
                                                       const su2double radius, unsigned long& maxErrorNodeLocal,
                                                       su2double& maxErrorLocal) {
  /*--- Array containing the local error ---*/
  su2double localError[3] = {0.0};

  /*--- Magnitude of the local maximum error ---*/
  maxErrorLocal = 0.0;

  /*--- Loop over non-selected boundary nodes ---*/
  for (auto iNode = 0ul; iNode < BoundNodes.size(); iNode++) {
    /*--- Compute nodal error ---*/
    GetNodalError(geometry, config, type, radius, iNode, localError);

    /*--- Setting error ---*/
    BoundNodes[iNode]->SetError(localError, nDim);

    /*--- Compute error magnitude and update local maximum error if necessary ---*/
    su2double errorMagnitude = GeometryToolbox::Norm(nDim, localError);
    if (errorMagnitude > maxErrorLocal) {
      maxErrorLocal = errorMagnitude;
      maxErrorNodeLocal = iNode;
    }
  }
}

void CRadialBasisFunctionInterpolation::GetNodalError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type,
                                                      const su2double radius, unsigned long iNode,
                                                      su2double* localError) {
  /*--- If requested (no by default) impose the surface deflections in increments ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  /*--- If node is part of a moving boundary then the error is defined as the difference
           between the found and prescribed displacements. Thus, here the displacement is substracted from the error
     ---*/
  if (config->GetMarker_All_Moving(BoundNodes[iNode]->GetMarker())) {
    auto displacement = geometry->vertex[BoundNodes[iNode]->GetMarker()][BoundNodes[iNode]->GetVertex()]->GetVarCoord();

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      localError[iDim] = -displacement[iDim] * VarIncrement;
    }
  } else {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      localError[iDim] = 0;
    }
  }

  /*--- Resulting displacement from the RBF interpolation is added to the error ---*/

  /*--- Finding contribution of each control node ---*/
  for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++) {
    /*--- Distance between non-selected boundary node and control node ---*/
    auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[jNode * nDim],
                                          geometry->nodes->GetCoord(BoundNodes[iNode]->GetIndex()));

    /*--- Evaluation of Radial Basis Function ---*/
    auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

    /*--- Add contribution to error ---*/
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      localError[iDim] += rbf * InterpCoeff[jNode * nDim + iDim];
    }
  }
}

void CRadialBasisFunctionInterpolation::SetCorrection(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type,
                                                      const vector<unsigned long>& internalNodes) {
  /*--- The non-selected control nodes still have a nonzero error once the maximum error falls below the data reduction
     tolerance. This error is applied as correction and interpolated into the volumetric mesh for internal nodes that
     fall within the correction radius. To evaluate whether an internal node falls within the correction radius an AD
     tree is constructed of the boundary nodes, making it possible to determine the distance to the nearest boundary
     node. ---*/

  /*--- Construction of the AD tree consisting of the non-selected boundary nodes ---*/

  /*--- Number of non-selected boundary nodes ---*/
  const unsigned long nVertexBound = BoundNodes.size();

  /*--- Vector storing the coordinates of the boundary nodes ---*/
  vector<su2double> Coord_bound(nDim * nVertexBound);

  /*--- Vector storing the IDs of the boundary nodes ---*/
  vector<unsigned long> PointIDs(nVertexBound);

  /*--- Correction Radius, equal to maximum error times a prescribed constant ---*/
  const su2double CorrectionRadius = config->GetRBFParam().GreedyCorrectionFactor * MaxErrorGlobal;

  /*--- Storing boundary node information ---*/
  unsigned long i = 0;
  unsigned long j = 0;
  for (auto iVertex = 0ul; iVertex < nVertexBound; iVertex++) {
    auto iNode = BoundNodes[iVertex]->GetIndex();
    PointIDs[i++] = iVertex;
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      Coord_bound[j++] = geometry->nodes->GetCoord(iNode, iDim);
    }
  }

  /*--- Construction of AD tree ---*/
  CADTPointsOnlyClass BoundADT(nDim, nVertexBound, Coord_bound.data(), PointIDs.data(), true);

  /*--- ID of nearest boundary node ---*/
  unsigned long pointID;
  /*--- Distance to nearest boundary node ---*/
  su2double dist;
  /*--- rank of nearest boundary node ---*/
  int rankID;

  /*--- Interpolation of the correction to the internal nodes that fall within the correction radius ---*/
  for (auto iNode = 0ul; iNode < internalNodes.size(); iNode++) {
    /*--- Find nearest node ---*/
    BoundADT.DetermineNearestNode(geometry->nodes->GetCoord(internalNodes[iNode]), dist, pointID, rankID);

    /*--- Get error of nearest node ---*/
    auto err = BoundNodes[pointID]->GetError();

    /*--- evaluate RBF ---*/
    auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, CorrectionRadius, dist));

    /*--- Apply correction to the internal node ---*/
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, -rbf * err[iDim]);
    }
  }

  /*--- Applying the correction to the non-selected boundary nodes ---*/
  for (auto iNode = 0ul; iNode < BoundNodes.size(); iNode++) {
    auto err = BoundNodes[iNode]->GetError();
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      geometry->nodes->AddCoord(BoundNodes[iNode]->GetIndex(), iDim, -err[iDim]);
    }
  }
}

void CRadialBasisFunctionInterpolation::AddControlNode(unsigned long maxErrorNode) {
  /*--- Addition of node to the reduced set of control nodes ---*/
  ReducedControlNodes.push_back(BoundNodes[maxErrorNode]);

  /*--- Removal of node among the non-selected boundary nodes ---*/
  BoundNodes.erase(BoundNodes.begin() + maxErrorNode);
}

void CRadialBasisFunctionInterpolation::ComputeNCtrlNodesGlobal() {
  /*--- Determining the global number of control nodes ---*/

  /*--- Local number of control nodes ---*/
  auto local_nControlNodes = ControlNodes->size();

  /*--- Summation of local number of control nodes ---*/
  SU2_MPI::Allreduce(&local_nControlNodes, &nCtrlNodesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
}
