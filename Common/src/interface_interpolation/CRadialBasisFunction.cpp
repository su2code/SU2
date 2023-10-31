/*!
 * \file CRadialBasisFunction.cpp
 * \brief Implementation of RBF interpolation.
 * \author Joel Ho, P. Gomes
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

#include "../../include/interface_interpolation/CRadialBasisFunction.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/CSymmetricMatrix.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

#if defined(HAVE_MKL)
#include "mkl.h"
#ifndef HAVE_LAPACK
#define HAVE_LAPACK
#endif
#elif defined(HAVE_LAPACK)
// dgemm(opA, opB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
extern "C" void dgemm_(const char*, const char*, const int*, const int*, const int*, const passivedouble*,
                       const passivedouble*, const int*, const passivedouble*, const int*, const passivedouble*,
                       passivedouble*, const int*);
#define DGEMM dgemm_
#endif

CRadialBasisFunction::CRadialBasisFunction(CGeometry**** geometry_container, const CConfig* const* config,
                                           unsigned int iZone, unsigned int jZone)
    : CInterpolator(geometry_container, config, iZone, jZone) {
  SetTransferCoeff(config);
}

void CRadialBasisFunction::PrintStatistics() const {
  if (rank != MASTER_NODE) return;
  cout.precision(3);
  cout << "  Min/avg/max number of RBF donors per target point: " << MinDonors << "/" << AvgDonors << "/" << MaxDonors
       << "\n"
       << "  Avg/max correction factor after pruning: " << AvgCorrection << "/" << MaxCorrection;
  if (MaxCorrection < 1.1 || AvgCorrection < 1.02)
    cout << " (ok)\n";
  else if (MaxCorrection < 2.0 && AvgCorrection < 1.05)
    cout << " (warning)\n";
  else
    cout << " <<< WARNING >>>\n";
  cout << "  Interpolation matrix is " << Density << "% dense." << endl;
  cout.unsetf(ios::floatfield);
}

su2double CRadialBasisFunction::Get_RadialBasisValue(RADIAL_BASIS type, const su2double radius, const su2double dist) {
  su2double rbf = dist / radius;

  switch (type) {
    case RADIAL_BASIS::WENDLAND_C2:
      if (rbf < 1)
        rbf = pow(pow((1 - rbf), 2), 2) * (4 * rbf + 1);  // double use of pow(x,2) for optimization
      else
        rbf = 0.0;
      break;

    case RADIAL_BASIS::GAUSSIAN:
      rbf = exp(-rbf * rbf);
      break;

    case RADIAL_BASIS::THIN_PLATE_SPLINE:
      if (rbf < numeric_limits<float>::min())
        rbf = 0.0;
      else
        rbf *= rbf * log(rbf);
      break;

    case RADIAL_BASIS::MULTI_QUADRIC:
    case RADIAL_BASIS::INV_MULTI_QUADRIC:
      rbf = sqrt(1.0 + rbf * rbf);
      if (type == RADIAL_BASIS::INV_MULTI_QUADRIC) rbf = 1.0 / rbf;
      break;
  }

  return rbf;
}

void CRadialBasisFunction::SetTransferCoeff(const CConfig* const* config) {
  /*--- RBF options. ---*/
  const auto kindRBF = config[donorZone]->GetKindRadialBasisFunction();
  const bool usePolynomial = config[donorZone]->GetRadialBasisFunctionPolynomialOption();
  const su2double paramRBF = config[donorZone]->GetRadialBasisFunctionParameter();
  const su2double pruneTol = config[donorZone]->GetRadialBasisFunctionPruneTol();

  const auto nMarkerInt = config[donorZone]->GetMarker_n_ZoneInterface() / 2;
  const int nDim = donor_geometry->GetnDim();

  const int nProcessor = size;
  Buffer_Receive_nVertex_Donor = new unsigned long[nProcessor];

  targetVertices.resize(config[targetZone]->GetnMarker_All());

  /*--- Process interface patches in parallel, fetch all donor point coordinates,
   *    then distribute interpolation matrix computation over ranks and threads.
   *    To avoid repeating calls to Collect_VertexInfo we also save the global
   *    indices of the donor points and the mpi rank index that owns them. ---*/

  vector<su2activematrix> donorCoordinates(nMarkerInt);
  vector<vector<long> > donorGlobalPoint(nMarkerInt);
  vector<vector<int> > donorProcessor(nMarkerInt);
  vector<int> assignedProcessor(nMarkerInt, -1);
  vector<unsigned long> totalWork(nProcessor, 0);

  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; ++iMarkerInt) {
    /*--- On the donor side: find the tag of the boundary sharing the interface. ---*/
    const auto markDonor = config[donorZone]->FindInterfaceMarker(iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface. ---*/
    const auto markTarget = config[targetZone]->FindInterfaceMarker(iMarkerInt);

    /*--- If the zone does not contain the interface continue to the next pair of markers. ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex(markDonor);

    /*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor. ---*/
    Determine_ArraySize(markDonor, markTarget, nVertexDonor, nDim);

    /*--- Compute total number of donor vertices. ---*/
    const auto nGlobalVertexDonor =
        accumulate(Buffer_Receive_nVertex_Donor, Buffer_Receive_nVertex_Donor + nProcessor, 0ul);

    /*--- Gather coordinates and global point indices. ---*/
    Buffer_Send_Coord.resize(MaxLocalVertex_Donor, nDim);
    Buffer_Send_GlobalPoint.resize(MaxLocalVertex_Donor);
    Buffer_Receive_Coord.resize(nProcessor * MaxLocalVertex_Donor, nDim);
    Buffer_Receive_GlobalPoint.resize(nProcessor * MaxLocalVertex_Donor);

    Collect_VertexInfo(markDonor, markTarget, nVertexDonor, nDim);

    /*--- Compresses the gathered donor point information to simplify computations. ---*/
    auto& donorCoord = donorCoordinates[iMarkerInt];
    auto& donorPoint = donorGlobalPoint[iMarkerInt];
    auto& donorProc = donorProcessor[iMarkerInt];
    donorCoord.resize(nGlobalVertexDonor, nDim);
    donorPoint.resize(nGlobalVertexDonor);
    donorProc.resize(nGlobalVertexDonor);

    auto iCount = 0ul;
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor) {
      auto offset = iProcessor * MaxLocalVertex_Donor;
      for (auto iVertex = 0ul; iVertex < Buffer_Receive_nVertex_Donor[iProcessor]; ++iVertex) {
        for (int iDim = 0; iDim < nDim; ++iDim) donorCoord(iCount, iDim) = Buffer_Receive_Coord(offset + iVertex, iDim);
        donorPoint[iCount] = Buffer_Receive_GlobalPoint[offset + iVertex];
        donorProc[iCount] = iProcessor;
        ++iCount;
      }
    }
    assert((iCount == nGlobalVertexDonor) && "Global donor point count mismatch.");

    /*--- Give an MPI-independent order to the points (required due to high condition
     *    number of the RBF matrix, avoids diff results with diff number of ranks. ---*/
    vector<int> order(nGlobalVertexDonor);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&donorPoint](int i, int j) { return donorPoint[i] < donorPoint[j]; });

    for (int i = 0; i < int(nGlobalVertexDonor); ++i) {
      int j = order[i];
      while (j < i) j = order[j];
      if (i == j) continue;
      swap(donorProc[i], donorProc[j]);
      swap(donorPoint[i], donorPoint[j]);
      for (int iDim = 0; iDim < nDim; ++iDim) swap(donorCoord(i, iDim), donorCoord(j, iDim));
    }

    /*--- Static work scheduling over ranks based on which one has less work currently. ---*/
    int iProcessor = 0;
    for (int i = 1; i < nProcessor; ++i)
      if (totalWork[i] < totalWork[iProcessor]) iProcessor = i;

    totalWork[iProcessor] += pow(nGlobalVertexDonor, 3);  // based on matrix inversion.

    assignedProcessor[iMarkerInt] = iProcessor;
  }
  delete[] Buffer_Receive_nVertex_Donor;

  /*--- Compute the interpolation matrices for each patch of coordinates
   *    assigned to the rank. Subdivide work further by threads. ---*/
  vector<int> nPolynomialVec(nMarkerInt, -1);
  vector<vector<int> > keepPolynomialRowVec(nMarkerInt, vector<int>(nDim, 1));
  vector<su2passivematrix> CinvTrucVec(nMarkerInt);

  SU2_OMP_PARALLEL_(for schedule(dynamic,1))
  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; ++iMarkerInt) {
    if (rank == assignedProcessor[iMarkerInt]) {
      ComputeGeneratorMatrix(kindRBF, usePolynomial, paramRBF, donorCoordinates[iMarkerInt], nPolynomialVec[iMarkerInt],
                             keepPolynomialRowVec[iMarkerInt], CinvTrucVec[iMarkerInt]);
    }
  }
  END_SU2_OMP_PARALLEL

  /*--- Final loop over interface markers to compute the interpolation coefficients. ---*/

  /*--- Initialize variables for interpolation statistics. ---*/
  unsigned long totalTargetPoints = 0, totalDonorPoints = 0, denseSize = 0;
  MinDonors = 1 << 30;
  MaxDonors = 0;
  MaxCorrection = 0.0;
  AvgCorrection = 0.0;

  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; iMarkerInt++) {
    /*--- Identify the rank that computed the interpolation matrix for this marker. ---*/
    const int iProcessor = assignedProcessor[iMarkerInt];
    /*--- If no processor was assigned to work, the zone does not contain the interface. ---*/
    if (iProcessor < 0) continue;

    /*--- Setup target information. ---*/
    const auto markTarget = config[targetZone]->FindInterfaceMarker(iMarkerInt);
    unsigned long nVertexTarget = 0;
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*--- Set references to donor information. ---*/
    auto& donorCoord = donorCoordinates[iMarkerInt];
    auto& donorPoint = donorGlobalPoint[iMarkerInt];
    auto& donorProc = donorProcessor[iMarkerInt];

    auto& C_inv_trunc = CinvTrucVec[iMarkerInt];
    auto& nPolynomial = nPolynomialVec[iMarkerInt];
    auto& keepPolynomialRow = keepPolynomialRowVec[iMarkerInt];

    const auto nGlobalVertexDonor = donorCoord.rows();

#ifdef HAVE_MPI
    /*--- For simplicity, broadcast small information about the interpolation matrix. ---*/
    SU2_MPI::Bcast(&nPolynomial, 1, MPI_INT, iProcessor, SU2_MPI::GetComm());
    SU2_MPI::Bcast(keepPolynomialRow.data(), nDim, MPI_INT, iProcessor, SU2_MPI::GetComm());

    /*--- Send C_inv_trunc only to the ranks that need it (those with target points),
     *    partial broadcast. MPI wrapper not used due to passive double. ---*/
    vector<unsigned long> allNumVertex(nProcessor);
    SU2_MPI::Allgather(&nVertexTarget, 1, MPI_UNSIGNED_LONG, allNumVertex.data(), 1, MPI_UNSIGNED_LONG,
                       SU2_MPI::GetComm());

    if (rank == iProcessor) {
      for (int jProcessor = 0; jProcessor < nProcessor; ++jProcessor)
        if ((jProcessor != iProcessor) && (allNumVertex[jProcessor] != 0))
          MPI_Send(C_inv_trunc.data(), C_inv_trunc.size(), MPI_DOUBLE, jProcessor, 0, SU2_MPI::GetComm());
    } else if (nVertexTarget != 0) {
      C_inv_trunc.resize(1 + nPolynomial + nGlobalVertexDonor, nGlobalVertexDonor);
      MPI_Recv(C_inv_trunc.data(), C_inv_trunc.size(), MPI_DOUBLE, iProcessor, 0, SU2_MPI::GetComm(),
               MPI_STATUS_IGNORE);
    }
#endif

    /*--- Compute interpolation matrix (H). This is a large matrix-matrix product with
     *    the generator matrix (C_inv_trunc) on the right. We avoid instantiation
     *    of the entire function matrix (A) and of the result (H), but work
     *    on a slab (set of rows) of A/H to amortize accesses to C_inv_trunc. ---*/

    /*--- Fetch target vertex coordinates. ---*/

    if (nVertexTarget) targetVertices[markTarget].resize(nVertexTarget);
    vector<const su2double*> targetCoord(nVertexTarget);

    for (auto iVertexTarget = 0ul; iVertexTarget < nVertexTarget; ++iVertexTarget) {
      const auto pointTarget = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();
      targetCoord[iVertexTarget] = target_geometry->nodes->GetCoord(pointTarget);
    }
    totalTargetPoints += nVertexTarget;
    denseSize += nVertexTarget * nGlobalVertexDonor;

    /*--- Distribute target slabs over the threads in the rank for processing. ---*/

    SU2_OMP_PARALLEL
    if (nVertexTarget > 0) {
      constexpr unsigned long targetSlabSize = 32;

      su2passivematrix funcMat(targetSlabSize, 1 + nPolynomial + nGlobalVertexDonor);
      su2passivematrix interpMat(targetSlabSize, nGlobalVertexDonor);

      /*--- Thread-local variables for statistics. ---*/
      unsigned long minDonors = 1 << 30, maxDonors = 0, totalDonors = 0;
      passivedouble sumCorr = 0.0, maxCorr = 0.0;

      SU2_OMP_FOR_DYN(1)
      for (auto iVertexTarget = 0ul; iVertexTarget < nVertexTarget; iVertexTarget += targetSlabSize) {
        const auto iLastVertex = min(nVertexTarget, iVertexTarget + targetSlabSize);
        const auto slabSize = iLastVertex - iVertexTarget;

        /*--- Prepare matrix of functions A (the targets to donors matrix). ---*/

        /*--- Polynominal part: ---*/
        if (usePolynomial) {
          /*--- Constant term. ---*/
          for (auto k = 0ul; k < slabSize; ++k) funcMat(k, 0) = 1.0;

          /*--- Linear terms. ---*/
          for (int iDim = 0, idx = 1; iDim < nDim; ++iDim) {
            /*--- Of which one may have been excluded. ---*/
            if (!keepPolynomialRow[iDim]) continue;
            for (auto k = 0ul; k < slabSize; ++k)
              funcMat(k, idx) = SU2_TYPE::GetValue(targetCoord[iVertexTarget + k][iDim]);
            idx += 1;
          }
        }
        /*--- RBF terms: ---*/
        for (auto iVertexDonor = 0ul; iVertexDonor < nGlobalVertexDonor; ++iVertexDonor) {
          for (auto k = 0ul; k < slabSize; ++k) {
            auto dist = GeometryToolbox::Distance(nDim, targetCoord[iVertexTarget + k], donorCoord[iVertexDonor]);
            auto rbf = Get_RadialBasisValue(kindRBF, paramRBF, dist);
            funcMat(k, 1 + nPolynomial + iVertexDonor) = SU2_TYPE::GetValue(rbf);
          }
        }

        /*--- Compute slab of the interpolation matrix. ---*/
#ifdef HAVE_LAPACK
        /*--- interpMat = funcMat * C_inv_trunc, but order of gemm arguments
         *    is swapped due to row-major storage of su2passivematrix. ---*/
        const char op = 'N';
        const int M = interpMat.cols(), N = slabSize, K = funcMat.cols();
        // lda = C_inv_trunc.cols() = M; ldb = funcMat.cols() = K; ldc = interpMat.cols() = M;
        const passivedouble alpha = 1.0, beta = 0.0;
        DGEMM(&op, &op, &M, &N, &K, &alpha, C_inv_trunc[0], &M, funcMat[0], &K, &beta, interpMat[0], &M);
#else
        /*--- Naive product, loop order considers short-wide
         *    nature of funcMat and interpMat. ---*/
        interpMat = 0.0;
        for (auto k = 0ul; k < funcMat.cols(); ++k)
          for (auto i = 0ul; i < slabSize; ++i)
            for (auto j = 0ul; j < interpMat.cols(); ++j) interpMat(i, j) += funcMat(i, k) * C_inv_trunc(k, j);
#endif
        /*--- Set interpolation coefficients. ---*/

        for (auto k = 0ul; k < slabSize; ++k) {
          auto& targetVertex = targetVertices[markTarget][iVertexTarget + k];

          /*--- Prune small coefficients. ---*/
          auto info = PruneSmallCoefficients(SU2_TYPE::GetValue(pruneTol), interpMat.cols(), interpMat[k]);
          auto nnz = info.first;
          totalDonors += nnz;
          minDonors = min(minDonors, nnz);
          maxDonors = max(maxDonors, nnz);
          auto corr = fabs(info.second - 1.0);  // far from 1 either way is bad;
          sumCorr += corr;
          maxCorr = max(maxCorr, corr);

          /*--- Allocate and set donor information for this target point. ---*/
          targetVertex.resize(nnz);

          for (unsigned long iVertex = 0, iSet = 0; iVertex < nGlobalVertexDonor; ++iVertex) {
            auto coeff = interpMat(k, iVertex);
            if (fabs(coeff) > 0.0) {
              targetVertex.processor[iSet] = donorProc[iVertex];
              targetVertex.globalPoint[iSet] = donorPoint[iVertex];
              targetVertex.coefficient[iSet] = coeff;
              ++iSet;
            }
          }
        }
      }  // end target vertex loop
      END_SU2_OMP_FOR
      SU2_OMP_CRITICAL {
        totalDonorPoints += totalDonors;
        MinDonors = min(MinDonors, minDonors);
        MaxDonors = max(MaxDonors, maxDonors);
        AvgCorrection += sumCorr;
        MaxCorrection = max(MaxCorrection, maxCorr);
      }
      END_SU2_OMP_CRITICAL
    }
    END_SU2_OMP_PARALLEL

    /*--- Free global data that will no longer be used. ---*/
    donorCoord.resize(0, 0);
    vector<long>().swap(donorPoint);
    vector<int>().swap(donorProc);
    C_inv_trunc.resize(0, 0);

  }  // end loop over interface markers

  /*--- Final reduction of interpolation statistics and basic sanity checks. ---*/
  auto Reduce = [](SU2_MPI::Op op, unsigned long& val) {
    auto tmp = val;
    SU2_MPI::Allreduce(&tmp, &val, 1, MPI_UNSIGNED_LONG, op, SU2_MPI::GetComm());
  };
  Reduce(MPI_SUM, totalTargetPoints);
  Reduce(MPI_SUM, totalDonorPoints);
  Reduce(MPI_SUM, denseSize);
  Reduce(MPI_MIN, MinDonors);
  Reduce(MPI_MAX, MaxDonors);
#ifdef HAVE_MPI
  passivedouble tmp1 = AvgCorrection, tmp2 = MaxCorrection;
  MPI_Allreduce(&tmp1, &AvgCorrection, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  MPI_Allreduce(&tmp2, &MaxCorrection, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
#endif
  if (totalTargetPoints == 0) SU2_MPI::Error("Somehow there are no target interpolation points.", CURRENT_FUNCTION);

  if (MinDonors == 0)
    SU2_MPI::Error(
        "One or more target points have no donors, either:\n"
        " - The interface surfaces are not in contact.\n"
        " - The RBF radius is too small.\n"
        " - The pruning tolerance is too aggressive.",
        CURRENT_FUNCTION);

  MaxCorrection += 1.0;  // put back the reference "1"
  AvgCorrection = AvgCorrection / totalTargetPoints + 1.0;
  AvgDonors = totalDonorPoints / totalTargetPoints;
  Density = totalDonorPoints / (0.01 * denseSize);
}

void CRadialBasisFunction::ComputeGeneratorMatrix(RADIAL_BASIS type, bool usePolynomial, su2double radius,
                                                  const su2activematrix& coords, int& nPolynomial,
                                                  vector<int>& keepPolynomialRow, su2passivematrix& C_inv_trunc) {
  const su2double interfaceCoordTol = 1e6 * numeric_limits<passivedouble>::epsilon();

  const int nVertexDonor = coords.rows();
  const int nDim = coords.cols();

  /*--- Populate interpolation kernel. ---*/
  CSymmetricMatrix global_M(nVertexDonor);

  for (int iVertex = 0; iVertex < nVertexDonor; ++iVertex)
    for (int jVertex = iVertex; jVertex < nVertexDonor; ++jVertex)
      global_M(iVertex, jVertex) = SU2_TYPE::GetValue(
          Get_RadialBasisValue(type, radius, GeometryToolbox::Distance(nDim, coords[iVertex], coords[jVertex])));

  /*--- Invert M matrix (operation is in-place). ---*/
  const bool kernelIsSPD = (type == RADIAL_BASIS::WENDLAND_C2) || (type == RADIAL_BASIS::GAUSSIAN) ||
                           (type == RADIAL_BASIS::INV_MULTI_QUADRIC);
  global_M.Invert(kernelIsSPD);

  /*--- Compute C_inv_trunc. ---*/
  if (usePolynomial) {
    /*--- Fill P matrix (P for points, with an extra top row of ones). ---*/
    su2passivematrix P(1 + nDim, nVertexDonor);

    for (int iVertex = 0; iVertex < nVertexDonor; iVertex++) {
      P(0, iVertex) = 1.0;
      for (int iDim = 0; iDim < nDim; ++iDim) P(1 + iDim, iVertex) = SU2_TYPE::GetValue(coords(iVertex, iDim));
    }

    /*--- Check if points lie on a plane and remove one coordinate from P if so. ---*/
    nPolynomial = CheckPolynomialTerms(interfaceCoordTol, keepPolynomialRow, P);

    /*--- Compute Q = P * M^-1 ---*/
    su2passivematrix Q;
    global_M.MatMatMult('R', P, Q);

    /*--- Compute Mp = (Q * P^T)^-1 ---*/
    CSymmetricMatrix Mp(nPolynomial + 1);

    for (int i = 0; i <= nPolynomial; ++i)
      for (int j = i; j <= nPolynomial; ++j) {
        Mp(i, j) = 0.0;
        for (int k = 0; k < nVertexDonor; ++k) Mp(i, j) += Q(i, k) * P(j, k);
      }
    Mp.Invert(false);

    /*--- Compute M_p * Q, the top part of C_inv_trunc. ---*/
    su2passivematrix C_inv_top;
    Mp.MatMatMult('L', Q, C_inv_top);

    /*--- Compute M^-1 - Q^T * C_inv_top, the bottom (main) part of C_inv_trunc. ---*/
    su2passivematrix C_inv_bot = global_M.StealData();
#ifdef HAVE_LAPACK
    /*--- Order of gemm arguments swapped due to row-major storage. ---*/
    const char opa = 'N', opb = 'T';
    const int M = nVertexDonor, N = nVertexDonor, K = nPolynomial + 1;
    // lda = C_inv_top.cols() = M; ldb = Q.cols() = M; ldc = C_inv_bot.cols() = M;
    const passivedouble alpha = -1.0, beta = 1.0;
    DGEMM(&opa, &opb, &M, &N, &K, &alpha, C_inv_top[0], &M, Q[0], &M, &beta, C_inv_bot[0], &M);
#else  // naive product
    for (int i = 0; i < nVertexDonor; ++i)
      for (int j = 0; j < nVertexDonor; ++j)
        for (int k = 0; k <= nPolynomial; ++k) C_inv_bot(i, j) -= Q(k, i) * C_inv_top(k, j);
#endif

    /*--- Merge top and bottom of C_inv_trunc. More intrusive memory
     *    management, or separate handling of top and bottom, would
     *    avoid these copies (and associated temporary vars). ---*/
    C_inv_trunc.resize(1 + nPolynomial + nVertexDonor, nVertexDonor);
    memcpy(C_inv_trunc[0], C_inv_top.data(), C_inv_top.size() * sizeof(passivedouble));
    memcpy(C_inv_trunc[1 + nPolynomial], C_inv_bot.data(), C_inv_bot.size() * sizeof(passivedouble));
  } else {
    /*--- No polynomial term used in the interpolation, C_inv_trunc = M^-1. ---*/

    C_inv_trunc = global_M.StealData();

  }  // end usePolynomial
}

int CRadialBasisFunction::CheckPolynomialTerms(su2double max_diff_tol, vector<int>& keep_row, su2passivematrix& P) {
  const int m = P.rows();
  const int n = P.cols();

  /*--- The first row of P is all ones and we do not care about it for this analysis. ---*/
  const int n_rows = m - 1;
  keep_row.resize(n_rows);

  /*--- By default assume points are not on a plane (all rows kept). ---*/
  int n_polynomial = n_rows;
  keep_row.resize(n_rows);
  for (int i = 0; i < n_rows; ++i) keep_row[i] = 1;

  /*--- Fit a plane through the points in P. ---*/

  /*--- Compute P times its transpose and invert. ---*/
  CSymmetricMatrix PPT(n_rows);

  for (int i = 0; i < n_rows; ++i)
    for (int j = i; j < n_rows; ++j) {
      PPT(i, j) = 0.0;
      for (int k = 0; k < n; ++k) PPT(i, j) += P(i + 1, k) * P(j + 1, k);
    }
  PPT.Invert(true);

  /*--- RHS for the least squares fit (vector of ones times P). ---*/
  vector<passivedouble> rhs(n_rows, 0.0), coeff(n_rows);

  for (int i = 0; i < n_rows; ++i)
    for (int j = 0; j < n; ++j) rhs[i] += P(i + 1, j);

  /*--- Multiply the RHS by the inverse thus obtaining the coefficients. ---*/
  PPT.MatVecMult(rhs.begin(), coeff.begin());

  /*--- Determine the maximum deviation of the points from the fitted plane. ---*/
  passivedouble max_diff = 0.0;

  for (int j = 0; j < n; ++j) {
    passivedouble sum = 0.0;
    for (int i = 0; i < n_rows; ++i) sum += coeff[i] * P(i + 1, j);

    /*--- 1.0 is the arbitrary constant we are assuming when fitting
     the plane, i.e. the vector of ones used to generate the RHS. ---*/
    max_diff = max(abs(1.0 - sum), max_diff);
  }

  /*--- If points lie on plane remove row associated with the maximum coefficient. ---*/
  if (max_diff < max_diff_tol) {
    /*--- Find the max coeff and mark the corresponding row for removal. ---*/
    int remove_row = 0;
    for (int i = 1; i < n_rows; ++i)
      if (abs(coeff[i]) > abs(coeff[remove_row])) remove_row = i;

    /*--- Mark row as removed and adjust number of polynomial terms. ---*/
    n_polynomial = n_rows - 1;
    keep_row[remove_row] = 0;

    /*--- Truncated P by shifting rows "up". ---*/
    for (auto i = remove_row + 1; i < m - 1; ++i)
      for (int j = 0; j < n; ++j) P(i, j) = P(i + 1, j);
  }

  return n_polynomial;
}
