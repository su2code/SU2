/*!
 * \file computeGradientsLeastSquares.hpp
 * \brief Generic implementation of Least-Squares gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
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

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

namespace detail {

/*!
 * \brief Prepare Smatrix for 2D.
 * \ingroup FvmAlgos
 */
FORCEINLINE void computeSmatrix(su2double r11, su2double r12, su2double r13,
                                su2double r22, su2double r23, su2double r33,
                                su2double detR2, su2double Smatrix[][2]) {
  Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
  Smatrix[0][1] = -r11*r12/detR2;
  Smatrix[1][1] = r11*r11/detR2;
}

/*!
 * \brief Prepare Smatrix for 3D.
 * \ingroup FvmAlgos
 */
FORCEINLINE void computeSmatrix(su2double r11, su2double r12, su2double r13,
                                su2double r22, su2double r23, su2double r33,
                                su2double detR2, su2double Smatrix[][3]) {
  su2double z11 = r22*r33;
  su2double z12 =-r12*r33;
  su2double z13 = r12*r23-r13*r22;
  su2double z22 = r11*r33;
  su2double z23 =-r11*r23;
  su2double z33 = r11*r22;

  Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
  Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
  Smatrix[0][2] = (z13*z33)/detR2;
  Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
  Smatrix[1][2] = (z23*z33)/detR2;
  Smatrix[2][2] = (z33*z33)/detR2;
}

/*!
 * \brief Solve the least-squares problem for one point.
 * \ingroup FvmAlgos
 * \note See detail::computeGradientsLeastSquares for the
 *       purpose of template "nDim" and "periodic".
 */
template<size_t nDim, bool periodic, class GradientType, class RMatrixType>
FORCEINLINE void solveLeastSquares(size_t iPoint,
                                   size_t varBegin,
                                   size_t varEnd,
                                   const RMatrixType& Rmatrix,
                                   GradientType& gradient)
{
  const auto eps = pow(std::numeric_limits<passivedouble>::epsilon(),2);

  /*--- Entries of upper triangular matrix R. ---*/

  if (periodic) {
    AD::StartPreacc();
    AD::SetPreaccIn(Rmatrix(iPoint,0,0));
    AD::SetPreaccIn(Rmatrix(iPoint,0,1));
    AD::SetPreaccIn(Rmatrix(iPoint,1,1));
  }

  su2double r11 = Rmatrix(iPoint,0,0);
  su2double r12 = Rmatrix(iPoint,0,1);
  su2double r22 = Rmatrix(iPoint,1,1);
  su2double r13 = 0.0, r23 = 0.0, r33 = 1.0;

  r11 = sqrt(max(r11, eps));
  r12 /= r11;
  r22 = sqrt(max(r22 - r12*r12, eps));

  if (nDim == 3) {
    if (periodic) {
      AD::SetPreaccIn(Rmatrix(iPoint,0,2));
      AD::SetPreaccIn(Rmatrix(iPoint,1,2));
      AD::SetPreaccIn(Rmatrix(iPoint,2,1));
      AD::SetPreaccIn(Rmatrix(iPoint,2,2));
    }

    r13 = Rmatrix(iPoint,0,2);
    r33 = Rmatrix(iPoint,2,2);
    const auto r23_a = Rmatrix(iPoint,1,2);
    const auto r23_b = Rmatrix(iPoint,2,1);

    r13 /= r11;
    r23 = r23_a/r22 - r23_b*r12/(r11*r22);
    r33 = sqrt(max(r33 - r23*r23 - r13*r13, eps));
  }

  /*--- Compute determinant ---*/

  const su2double detR2 = pow(r11*r22*r33, 2);

  /*--- S matrix := inv(R)*traspose(inv(R)) ---*/

  su2double Smatrix[nDim][nDim] = {{0.0}};

  /*--- Detect singular matrix ---*/

  if (detR2 > eps) {
    computeSmatrix(r11, r12, r13, r22, r23, r33, detR2, Smatrix);
  }

  if (periodic) {
    /*--- Stop preacc here as gradient is in/out. ---*/
    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = iDim; jDim < nDim; ++jDim)
        AD::SetPreaccOut(Smatrix[iDim][jDim]);
    AD::EndPreacc();
  }

  /*--- Computation of the gradient: S*c ---*/

  for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
  {
    su2double Cvector[nDim] = {0.0};

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = 0; jDim < nDim; ++jDim)
        Cvector[iDim] += Smatrix[min(iDim,jDim)][max(iDim,jDim)] * gradient(iPoint, iVar, jDim);

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      gradient(iPoint, iVar, iDim) = Cvector[iDim];
  }

  if (!periodic) {
    /*--- Stop preacc here instead as gradient is only out. ---*/
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        AD::SetPreaccOut(gradient(iPoint, iVar, iDim));
    AD::EndPreacc();
  }
}

/*!
 * \brief Compute the gradient of a field using inverse-distance-weighted or
 *        unweighted Least-Squares approximation.
 * \ingroup FvmAlgos
 * \note See notes from computeGradientsGreenGauss.hpp.
 * \param[in] solver - Optional, solver associated with the field (used only for MPI).
 * \param[in] kindMpiComm - Type of MPI communication required.
 * \param[in] kindPeriodicComm - Type of periodic communication required.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] weighted - Use inverse-distance weights.
 * \param[in] config - Configuration of the problem, used to identify types of boundaries.
 * \param[in] field - Generic object implementing operator (iPoint, iVar).
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 * \param[out] Rmatrix - Generic object implementing operator (iPoint, iDim, iDim).
 */
template<size_t nDim, class FieldType, class GradientType, class RMatrixType>
void computeGradientsLeastSquares(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  const CConfig& config,
                                  bool weighted,
                                  const FieldType& field,
                                  size_t varBegin,
                                  size_t varEnd,
                                  GradientType& gradient,
                                  RMatrixType& Rmatrix)
{
  const bool periodic = (solver != nullptr) && (config.GetnMarker_Periodic() > 0);

  const size_t nPointDomain = geometry.GetnPointDomain();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  size_t chunkSize = computeStaticChunkSize(nPointDomain,
                     omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- First loop over non-halo points of the grid. ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
  {
    auto nodes = geometry.nodes;
    const auto coord_i = nodes->GetCoord(iPoint);

    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();
    AD::SetPreaccIn(coord_i, nDim);

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      AD::SetPreaccIn(field(iPoint,iVar));

    /*--- Clear gradient and Rmatrix. ---*/

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = 0.0;

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = 0; jDim < nDim; ++jDim)
        Rmatrix(iPoint, iDim, jDim) = 0.0;


    for (auto jPoint : nodes->GetPoints(iPoint))
    {
      const auto coord_j = geometry.nodes->GetCoord(jPoint);
      AD::SetPreaccIn(coord_j, nDim);


      /*--- Distance vector from iPoint to jPoint ---*/

      su2double dist_ij[nDim] = {0.0};
      GeometryToolbox::Distance(nDim, coord_j, coord_i, dist_ij);


      /*--- Compute inverse weight, default 1 (unweighted). ---*/

      su2double weight = 1.0;
      if(weighted) weight = GeometryToolbox::SquaredNorm(nDim, dist_ij);

      /*--- Summations for entries of upper triangular matrix R. ---*/

      if (weight > 0.0)
      {
        weight = 1.0 / weight;

        for (size_t iDim = 0; iDim < nDim; ++iDim)
          for (size_t jDim = iDim; jDim < nDim; ++jDim)
            Rmatrix(iPoint,iDim,jDim) += dist_ij[iDim]*dist_ij[jDim]*weight;

        if (nDim == 3)
          Rmatrix(iPoint,2,1) += dist_ij[0]*dist_ij[nDim-1]*weight;

        /*--- Entries of c:= transpose(A)*b ---*/

        for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        {
          AD::SetPreaccIn(field(jPoint,iVar));

          su2double delta_ij = weight * (field(jPoint,iVar) - field(iPoint,iVar));

          for (size_t iDim = 0; iDim < nDim; ++iDim)
            gradient(iPoint, iVar, iDim) += dist_ij[iDim] * delta_ij;
        }
      }
    }

    if (periodic)
    {
      /*--- A second loop is required after periodic comms, checkpoint the preacc. ---*/

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        for (size_t jDim = 0; jDim < nDim; ++jDim)
          AD::SetPreaccOut(Rmatrix(iPoint, iDim, jDim));

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        for (size_t iDim = 0; iDim < nDim; ++iDim)
          AD::SetPreaccOut(gradient(iPoint, iVar, iDim));

      AD::EndPreacc();
    }
    else {
      /*--- Periodic comms are not needed, solve the LS problem for iPoint. ---*/

      solveLeastSquares<nDim, false>(iPoint, varBegin, varEnd, Rmatrix, gradient);
    }
  }
  END_SU2_OMP_FOR

  /*--- Correct the gradient values across any periodic boundaries. ---*/

  if (periodic)
  {
    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
    {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    }

    /*--- Second loop over points of the grid to compute final gradient. ---*/

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
      solveLeastSquares<nDim, true>(iPoint, varBegin, varEnd, Rmatrix, gradient);
    END_SU2_OMP_FOR
  }

  /*--- If no solver was provided we do not communicate ---*/

  if (solver != nullptr)
  {
    /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

    solver->InitiateComms(&geometry, &config, kindMpiComm);
    solver->CompleteComms(&geometry, &config, kindMpiComm);
  }

}
} // end namespace

/*!
 * \brief Instantiations for 2D and 3D.
 * \ingroup FvmAlgos
 */
template<class FieldType, class GradientType, class RMatrixType>
void computeGradientsLeastSquares(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  const CConfig& config,
                                  bool weighted,
                                  const FieldType& field,
                                  size_t varBegin,
                                  size_t varEnd,
                                  GradientType& gradient,
                                  RMatrixType& Rmatrix) {
  switch (geometry.GetnDim()) {
  case 2:
    detail::computeGradientsLeastSquares<2>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                            weighted, field, varBegin, varEnd, gradient, Rmatrix);
    break;
  case 3:
    detail::computeGradientsLeastSquares<3>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                            weighted, field, varBegin, varEnd, gradient, Rmatrix);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}
