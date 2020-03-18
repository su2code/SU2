/*!
 * \file CRadialBasisFunction.cpp
 * \brief Implementation of RBF interpolation.
 * \author Joel Ho, P. Gomes
 * \version 7.0.2 "Blackbird"
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

#include "../../include/interface_interpolation/CRadialBasisFunction.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/CSymmetricMatrix.hpp"


CRadialBasisFunction::CRadialBasisFunction(CGeometry ****geometry_container, CConfig **config, unsigned int iZone,
                                           unsigned int jZone) : CInterpolator(geometry_container, config, iZone, jZone) {
  Set_TransferCoeff(config);
}

su2double CRadialBasisFunction::Get_RadialBasisValue(ENUM_RADIALBASIS type, const su2double radius, const su2double dist)
{
  su2double rbf = dist/radius;

  switch (type) {

    case WENDLAND_C2:
      if(rbf < 1) rbf = pow(pow((1-rbf),2),2)*(4*rbf+1); // double use of pow(x,2) for optimization
      else        rbf = 0.0;
      break;

    case GAUSSIAN:
      rbf = exp(-rbf*rbf);
      break;

    case THIN_PLATE_SPLINE:
      if(rbf < numeric_limits<float>::min()) rbf = 0.0;
      else rbf *= rbf*log(rbf);
      break;

    case MULTI_QUADRIC:
    case INV_MULTI_QUADRIC:
      rbf = sqrt(1.0+rbf*rbf);
      if(type == INV_MULTI_QUADRIC) rbf = 1.0/rbf;
      break;
  }

  return rbf;
}

void CRadialBasisFunction::Set_TransferCoeff(CConfig **config) {

  /*--- RBF options. ---*/
  const auto kindRBF = static_cast<ENUM_RADIALBASIS>(config[donorZone]->GetKindRadialBasisFunction());
  const bool usePolynomial = config[donorZone]->GetRadialBasisFunctionPolynomialOption();
  const su2double paramRBF = config[donorZone]->GetRadialBasisFunctionParameter();
  const su2double pruneTol = config[donorZone]->GetRadialBasisFunctionPruneTol();

  const auto nMarkerInt = config[donorZone]->GetMarker_n_ZoneInterface()/2;
  const int nDim = donor_geometry->GetnDim();

  const int nProcessor = size;
  Buffer_Send_nVertex_Donor = new unsigned long [1];
  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];

  /*--- Process interface patches in parallel, fetch all donor point coordinates,
   *    then distribute interpolation matrix computation over ranks and threads.
   *    To avoid repeating calls to Collect_VertexInfo we also save the global
   *    indices of the donor points and the mpi rank index that owns them. ---*/
  vector<su2activematrix> DonorCoordinates(nMarkerInt);
  vector<vector<long> > DonorGlobalPoint(nMarkerInt);
  vector<vector<int> > DonorProcessor(nMarkerInt);
  vector<int> AssignedProcessor(nMarkerInt,-1);
  vector<unsigned long> TotalWork(nProcessor,0);

  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; ++iMarkerInt) {

    /*--- On the donor side: find the tag of the boundary sharing the interface. ---*/
    const auto mark_donor = Find_InterfaceMarker(config[donorZone], iMarkerInt+1);

    /*--- On the target side: find the tag of the boundary sharing the interface. ---*/
    const auto mark_target = Find_InterfaceMarker(config[targetZone], iMarkerInt+1);

    /*--- If the zone does not contain the interface continue to the next pair of markers. ---*/
    if(!CheckInterfaceBoundary(mark_donor,mark_target)) continue;

    unsigned long nVertexDonor = 0;
    if(mark_donor != -1) nVertexDonor = donor_geometry->GetnVertex(mark_donor);

    /*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor. ---*/
    Determine_ArraySize(false, mark_donor, mark_target, nVertexDonor, nDim);

    /*--- Compute total number of donor vertices. ---*/
    auto nGlobalVertexDonor = accumulate(Buffer_Receive_nVertex_Donor,
                              Buffer_Receive_nVertex_Donor+nProcessor, 0ul);

    /*--- Gather coordinates and global point indices. ---*/
    Buffer_Send_Coord = new su2double [ MaxLocalVertex_Donor * nDim ];
    Buffer_Send_GlobalPoint = new long [ MaxLocalVertex_Donor ];
    Buffer_Receive_Coord = new su2double [ nProcessor * MaxLocalVertex_Donor * nDim ];
    Buffer_Receive_GlobalPoint = new long [ nProcessor * MaxLocalVertex_Donor ];

    Collect_VertexInfo(false, mark_donor, mark_target, nVertexDonor, nDim);

    /*--- Compresses the gathered donor point information to simplify computations. ---*/
    auto& DonorCoord = DonorCoordinates[iMarkerInt];
    auto& DonorPoint = DonorGlobalPoint[iMarkerInt];
    auto& DonorProc = DonorProcessor[iMarkerInt];
    DonorCoord.resize(nGlobalVertexDonor, nDim);
    DonorPoint.resize(nGlobalVertexDonor);
    DonorProc.resize(nGlobalVertexDonor);

    auto iCount = 0ul;
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor) {
      auto offset = iProcessor * MaxLocalVertex_Donor;
      for (auto iVertex = 0ul; iVertex < Buffer_Receive_nVertex_Donor[iProcessor]; ++iVertex) {
        for (int iDim = 0; iDim < nDim; ++iDim)
          DonorCoord(iCount,iDim) = Buffer_Receive_Coord[(offset+iVertex)*nDim + iDim];
        DonorPoint[iCount] = Buffer_Receive_GlobalPoint[offset+iVertex];
        DonorProc[iCount] = iProcessor;
        ++iCount;
      }
    }
    assert((iCount == nGlobalVertexDonor) && "Global donor point count mismatch.");

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_GlobalPoint;
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_GlobalPoint;

    /*--- Static work scheduling over ranks based on which one has less work currently. ---*/
    int iProcessor = 0;
    for (int i = 1; i < nProcessor; ++i)
      if (TotalWork[i] < TotalWork[iProcessor]) iProcessor = i;

    TotalWork[iProcessor] += pow(nGlobalVertexDonor,3); // based on matrix inversion.

    AssignedProcessor[iMarkerInt] = iProcessor;

  }

  /*--- Compute the interpolation matrices for each patch of coordinates
   *    assigned to the rank. Subdivide work further by threads. ---*/
  vector<int> nPolynomialVec(nMarkerInt,-1);
  vector<vector<int> > keepPolynomialRowVec(nMarkerInt, vector<int>(nDim,1));
  vector<su2passivematrix> CinvTrucVec(nMarkerInt);

  SU2_OMP_PARALLEL_(for schedule(dynamic,1))
  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; ++iMarkerInt) {
    if (rank == AssignedProcessor[iMarkerInt]) {
      ComputeGeneratorMatrix(kindRBF, usePolynomial, paramRBF,
                             DonorCoordinates[iMarkerInt], nPolynomialVec[iMarkerInt],
                             keepPolynomialRowVec[iMarkerInt], CinvTrucVec[iMarkerInt]);
    }
  }

  /*--- Final loop over interface markers to compute the interpolation coefficients. ---*/

  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; iMarkerInt++) {

    /*--- Identify the rank that computed the interpolation matrix for this marker. ---*/
    const int iProcessor = AssignedProcessor[iMarkerInt];
    /*--- If no processor was assigned to work, the zone does not contain the interface. ---*/
    if (iProcessor < 0) continue;

    /*--- Setup target information. ---*/
    const int mark_target = Find_InterfaceMarker(config[targetZone], iMarkerInt+1);
    unsigned long nVertexTarget = 0;
    if(mark_target != -1) nVertexTarget = target_geometry->GetnVertex(mark_target);

    /*--- Set references to donor information. ---*/
    auto& DonorCoord = DonorCoordinates[iMarkerInt];
    auto& DonorPoint = DonorGlobalPoint[iMarkerInt];
    auto& DonorProc = DonorProcessor[iMarkerInt];

    auto& C_inv_trunc = CinvTrucVec[iMarkerInt];
    auto& nPolynomial = nPolynomialVec[iMarkerInt];
    auto& keepPolynomialRow = keepPolynomialRowVec[iMarkerInt];

    const auto nGlobalVertexDonor = DonorCoord.rows();

#ifdef HAVE_MPI
    /*--- For simplicity, broadcast small information about the interpolation matrix. ---*/
    SU2_MPI::Bcast(&nPolynomial, 1, MPI_INT, iProcessor, MPI_COMM_WORLD);
    SU2_MPI::Bcast(keepPolynomialRow.data(), nDim, MPI_INT, iProcessor, MPI_COMM_WORLD);

    /*--- Send C_inv_trunc only to the ranks that need it (those with target points),
     *    partial broadcast. MPI wrapper not used due to passive double. ---*/
    vector<unsigned long> allNumVertex(nProcessor);
    SU2_MPI::Allgather(&nVertexTarget, 1, MPI_UNSIGNED_LONG,
      allNumVertex.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    if (rank == iProcessor) {
      for (int jProcessor = 0; jProcessor < nProcessor; ++jProcessor)
        if ((jProcessor != iProcessor) && (allNumVertex[jProcessor] != 0))
          MPI_Send(C_inv_trunc.data(), C_inv_trunc.size(),
                   MPI_DOUBLE, jProcessor, 0, MPI_COMM_WORLD);
    }
    else if (nVertexTarget != 0) {
      C_inv_trunc.resize(1+nPolynomial+nGlobalVertexDonor, nGlobalVertexDonor);
      MPI_Recv(C_inv_trunc.data(), C_inv_trunc.size(), MPI_DOUBLE,
               iProcessor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
#endif

    /*--- Compute H (interpolation) matrix, distributing
     *    target points over the threads in the rank. ---*/
    SU2_OMP_PARALLEL
    {
    su2passivevector coeff_vec(nGlobalVertexDonor);

    SU2_OMP_FOR_DYN(roundUpDiv(nVertexTarget, 2*omp_get_max_threads()))
    for (auto iVertexTarget = 0ul; iVertexTarget < nVertexTarget; iVertexTarget++) {

      auto target_vertex = target_geometry->vertex[mark_target][iVertexTarget];
      const auto point_target = target_vertex->GetNode();

      /*--- If not domain point move to next. ---*/
      if (!target_geometry->node[point_target]->GetDomain()) continue;

      const su2double* coord_i = target_geometry->node[point_target]->GetCoord();

      /*--- Multiply target vector by C_inv_trunc to obtain the interpolation coefficients.
       *    The target vector is not stored, we consume its entries immediately to avoid
       *    strided access to C_inv_trunc (as it is row-major). ---*/

      /*--- Polynominal part: ---*/
      if (usePolynomial) {
        /*--- Constant term. ---*/
        for (auto j = 0ul; j < nGlobalVertexDonor; ++j)
          coeff_vec(j) = C_inv_trunc(0,j);

        /*--- Linear terms. ---*/
        for (int iDim = 0, idx = 1; iDim < nDim; ++iDim) {
          /*--- Of which one may have been excluded. ---*/
          if (!keepPolynomialRow[iDim]) continue;
          for (auto j = 0ul; j < nGlobalVertexDonor; ++j)
            coeff_vec(j) += SU2_TYPE::GetValue(coord_i[iDim]) * C_inv_trunc(idx,j);
          idx += 1;
        }
      }
      else {
        /*--- Initialize vector to zero. ---*/
        for (auto j = 0ul; j < nGlobalVertexDonor; ++j) coeff_vec(j) = 0.0;
      }

      /*--- RBF terms: ---*/
      for (auto iVertexDonor = 0ul; iVertexDonor < nGlobalVertexDonor; ++iVertexDonor) {
        auto w_ij = SU2_TYPE::GetValue(Get_RadialBasisValue(kindRBF, paramRBF,
                    PointsDistance(nDim, coord_i, DonorCoord[iVertexDonor])));

        for (auto j = 0ul; j < nGlobalVertexDonor; ++j)
          coeff_vec(j) += w_ij * C_inv_trunc(1+nPolynomial+iVertexDonor, j);
      }

      /*--- Prune small coefficients. ---*/
      auto nnz = PruneSmallCoefficients(SU2_TYPE::GetValue(pruneTol), coeff_vec);

      /*--- Allocate and set donor information for this target point. ---*/
      target_vertex->SetnDonorPoints(nnz);
      target_vertex->Allocate_DonorInfo();

      for (unsigned long iVertex = 0, iSet = 0; iVertex < nGlobalVertexDonor; ++iVertex) {
        if (fabs(coeff_vec(iVertex)) > 0.0) {
          target_vertex->SetInterpDonorProcessor(iSet, DonorProc[iVertex]);
          target_vertex->SetInterpDonorPoint(iSet, DonorPoint[iVertex]);
          target_vertex->SetDonorCoeff(iSet, coeff_vec(iVertex));
          ++iSet;
        }
      }

    } // end target vertex loop
    } // end SU2_OMP_PARALLEL

    /*--- Delete global data that will no longer be used. ---*/
    DonorCoord.resize(0,0);
    vector<long>().swap(DonorPoint);
    vector<int>().swap(DonorProc);
    C_inv_trunc.resize(0,0);

  } // end loop over interface markers

  delete[] Buffer_Send_nVertex_Donor;
  delete[] Buffer_Receive_nVertex_Donor;

}

void CRadialBasisFunction::ComputeGeneratorMatrix(ENUM_RADIALBASIS type, bool usePolynomial,
                           su2double radius, const su2activematrix& coords, int& nPolynomial,
                           vector<int>& keepPolynomialRow, su2passivematrix& C_inv_trunc) {

  const su2double interfaceCoordTol = 1e6 * numeric_limits<passivedouble>::epsilon();

  const auto nVertexDonor = coords.rows();
  const int nDim = coords.cols();

  /*--- Populate interpolation kernel. ---*/
  CSymmetricMatrix global_M(nVertexDonor);

  for (auto iVertex = 0ul; iVertex < nVertexDonor; ++iVertex)
    for (auto jVertex = iVertex; jVertex < nVertexDonor; ++jVertex)
      global_M(iVertex, jVertex) = SU2_TYPE::GetValue(Get_RadialBasisValue(type, radius,
                                   PointsDistance(nDim, coords[iVertex], coords[jVertex])));

  /*--- Invert M matrix (operation is in-place). ---*/
  const bool kernelIsSPD = (type==WENDLAND_C2) || (type==GAUSSIAN) || (type==INV_MULTI_QUADRIC);
  global_M.Invert(kernelIsSPD);

  /*--- Compute C_inv_trunc. ---*/
  if (usePolynomial) {

    /*--- Fill P matrix (P for points, with an extra top row of ones). ---*/
    su2passivematrix P(1+nDim, nVertexDonor);

    for (auto iVertex = 0ul; iVertex < nVertexDonor; iVertex++) {
      P(0, iVertex) = 1.0;
      for (int iDim = 0; iDim < nDim; ++iDim)
        P(1+iDim, iVertex) = SU2_TYPE::GetValue(coords(iVertex, iDim));
    }

    /*--- Check if points lie on a plane and remove one coordinate from P if so. ---*/
    nPolynomial = CheckPolynomialTerms(interfaceCoordTol, keepPolynomialRow, P);

    /*--- Compute Mp = (P * M^-1 * P^T)^-1 ---*/
    CSymmetricMatrix Mp(nPolynomial+1);

    su2passivematrix tmp;
    global_M.MatMatMult('R', P, tmp); // tmp = P * M^-1

    for (int i = 0; i <= nPolynomial; ++i) // Mp = tmp * P
      for (int j = i; j <= nPolynomial; ++j) {
        Mp(i,j) = 0.0;
        for (auto k = 0ul; k < nVertexDonor; ++k) Mp(i,j) += tmp(i,k) * P(j,k);
      }
    Mp.Invert(false); // Mp = Mp^-1

    /*--- Compute M_p * P * M^-1, the top part of C_inv_trunc. ---*/
    Mp.MatMatMult('L', P, tmp);
    su2passivematrix C_inv_top;
    global_M.MatMatMult('R', tmp, C_inv_top);

    /*--- Compute tmp = (I - P^T * M_p * P * M^-1), part of the bottom part of
     C_inv_trunc. Note that most of the product is known from the top part. ---*/
    tmp.resize(nVertexDonor, nVertexDonor);

    for (auto i = 0ul; i < nVertexDonor; ++i) {
      for (auto j = 0ul; j < nVertexDonor; ++j) {
        tmp(i,j) = 0.0;
        for (int k = 0; k <= nPolynomial; ++k) tmp(i,j) -= P(k,i) * C_inv_top(k,j);
      }
      tmp(i,i) += 1.0; // identity part
    }

    /*--- Compute M^-1 * (I - P^T * M_p * P * M^-1), finalize bottom of C_inv_trunc. ---*/
    global_M.MatMatMult('L', tmp, C_inv_trunc);

    /*--- Merge top and bottom of C_inv_trunc. ---*/
    tmp = move(C_inv_trunc);
    C_inv_trunc.resize(1+nPolynomial+nVertexDonor, nVertexDonor);
    memcpy(C_inv_trunc[0], C_inv_top.data(), C_inv_top.size()*sizeof(passivedouble));
    memcpy(C_inv_trunc[1+nPolynomial], tmp.data(), tmp.size()*sizeof(passivedouble));
  }
  else {
    /*--- No polynomial term used in the interpolation, C_inv_trunc = M^-1. ---*/

    C_inv_trunc.resize(nVertexDonor, nVertexDonor);
    for (auto i = 0ul; i < nVertexDonor; ++i)
      for (auto j = 0ul; j < nVertexDonor; ++j)
        C_inv_trunc(i,j) = global_M(i,j);

  } // end usePolynomial

}

int CRadialBasisFunction::CheckPolynomialTerms(su2double max_diff_tol, vector<int>& keep_row,
                                               su2passivematrix &P) {
  const int m = P.rows();
  const int n = P.cols();

  /*--- The first row of P is all ones and we do not care about it for this analysis. ---*/
  const int n_rows = m-1;
  keep_row.resize(n_rows);

  /*--- By default assume points are not on a plane (all rows kept). ---*/
  int n_polynomial = n_rows;
  for (int i = 0; i < n_rows; ++i) keep_row[i] = 1;

  /*--- Fit a plane through the points in P. ---*/

  /*--- Compute P times its transpose and invert. ---*/
  CSymmetricMatrix PPT(n_rows);

  for (int i = 0; i < n_rows; ++i)
    for (int j = i; j < n_rows; ++j) {
      PPT(i,j) = 0.0;
      for (int k = 0; k < n; ++k) PPT(i,j) += P(i+1,k) * P(j+1,k);
    }
  PPT.Invert(true);

  /*--- RHS for the least squares fit (vector of ones times P). ---*/
  vector<passivedouble> coeff(n_rows,0.0);

  for (int i = 0; i < n_rows; ++i)
    for (int j = 0; j < n; ++j)
      coeff[i] += P(i+1,j);

  /*--- Multiply the RHS by the inverse thus obtaining the coefficients. ---*/
  PPT.MatVecMult(coeff.data());

  /*--- Determine the maximum deviation of the points from the fitted plane. ---*/
  passivedouble max_diff = 0.0;

  for (int j = 0; j < n; ++j)
  {
    passivedouble sum = 0.0;
    for (int i = 0; i < n_rows; ++i) sum += coeff[i] * P(i+1,j);

    /*--- 1.0 is the arbitrary constant we are assuming when fitting
     the plane, i.e. the vector of ones used to generate the RHS. ---*/
    max_diff = max(abs(1.0-sum), max_diff);
  }

  /*--- If points lie on plane remove row associated with the maximum coefficient. ---*/
  if (max_diff < max_diff_tol)
  {
    /*--- Find the max coeff and mark the corresponding row for removal. ---*/
    int remove_row = 0;
    for (int i = 1; i < n_rows; ++i)
      if (abs(coeff[i]) > abs(coeff[remove_row]))
        remove_row = i;

    /*--- Mark row as removed and adjust number of polynomial terms. ---*/
    n_polynomial = n_rows-1;
    keep_row[remove_row] = 0;

    /*--- Truncated P by shifting rows "up". ---*/
    for (auto i = remove_row+1; i < m-1; ++i)
      for (int j = 0; j < n; ++j)
        P(i,j) = P(i+1,j);
  }

  return n_polynomial;
}

int CRadialBasisFunction::PruneSmallCoefficients(passivedouble tolerance,
                                                 su2passivevector& coeffs) {

  /*--- Determine the pruning threshold. ---*/
  passivedouble thresh = 0.0;
  for (auto i = 0ul; i < coeffs.size(); ++i)
    thresh = max(thresh, fabs(coeffs(i)));
  thresh *= tolerance;

  /*--- Prune and count non-zeros. ---*/
  int numNonZeros = 0;
  passivedouble coeffSum = 0.0;
  for (auto i = 0ul; i < coeffs.size(); ++i) {
    if (fabs(coeffs(i)) > thresh) {
      coeffSum += coeffs(i);
      ++numNonZeros;
    }
    else coeffs(i) = 0.0;
  }

  /*--- Correct remaining coefficients, sum must be 1 for conservation. ---*/
  passivedouble correction = 1.0 / coeffSum;
  for (auto i = 0ul; i < coeffs.size(); ++i) coeffs(i) *= correction;

  return numNonZeros;
}
