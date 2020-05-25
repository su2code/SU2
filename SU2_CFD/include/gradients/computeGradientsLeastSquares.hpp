/*!
 * \file computeGradientsLeastSquares.hpp
 * \brief Generic implementation of Least-Squares gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
 * \version 7.0.1 "Blackbird"
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

#include "../../../Common/include/omp_structure.hpp"


/*!
 * \brief Compute the gradient of a field using inverse-distance-weighted or
 *        unweighted Least-Squares approximation.
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
template<class FieldType, class GradientType, class RMatrixType>
void computeGradientsLeastSquares(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  CConfig& config,
                                  bool weighted,
                                  const FieldType& field,
                                  size_t varBegin,
                                  size_t varEnd,
                                  GradientType& gradient,
                                  RMatrixType& Rmatrix)
{
  constexpr size_t MAXNDIM = 3;

  size_t nPointDomain = geometry.GetnPointDomain();
  size_t nDim = geometry.GetnDim();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  size_t chunkSize = computeStaticChunkSize(nPointDomain,
                     omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- Start OpenMP parallel section. ---*/

  SU2_OMP_PARALLEL
  {
    /*--- First loop over non-halo points of the grid. ---*/

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
    {
      auto node = geometry.node[iPoint];
      const su2double* coord_i = node->GetCoord();

      AD::StartPreacc();
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


      for (size_t iNeigh = 0; iNeigh < node->GetnPoint(); ++iNeigh)
      {
        size_t jPoint = node->GetPoint(iNeigh);

        const su2double* coord_j = geometry.node[jPoint]->GetCoord();
        AD::SetPreaccIn(coord_j, nDim);

        /*--- Distance vector from iPoint to jPoint ---*/

        su2double dist_ij[MAXNDIM] = {0.0};

        for (size_t iDim = 0; iDim < nDim; ++iDim)
          dist_ij[iDim] = coord_j[iDim] - coord_i[iDim];

        /*--- Compute inverse weight, default 1 (unweighted). ---*/

        su2double weight = 1.0;

        if (weighted)
        {
          weight = 0.0;
          for (size_t iDim = 0; iDim < nDim; ++iDim)
            weight += dist_ij[iDim] * dist_ij[iDim];
        }

        /*--- Sumations for entries of upper triangular matrix R. ---*/

        if (weight > 0.0)
        {
          weight = 1.0 / weight;

          Rmatrix(iPoint,0,0) += dist_ij[0]*dist_ij[0]*weight;
          Rmatrix(iPoint,0,1) += dist_ij[0]*dist_ij[1]*weight;
          Rmatrix(iPoint,1,1) += dist_ij[1]*dist_ij[1]*weight;

          if (nDim == 3)
          {
            Rmatrix(iPoint,0,2) += dist_ij[0]*dist_ij[2]*weight;
            Rmatrix(iPoint,1,2) += dist_ij[1]*dist_ij[2]*weight;
            Rmatrix(iPoint,2,1) += dist_ij[0]*dist_ij[2]*weight;
            Rmatrix(iPoint,2,2) += dist_ij[2]*dist_ij[2]*weight;
          }

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

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        for (size_t jDim = 0; jDim < nDim; ++jDim)
          AD::SetPreaccOut(Rmatrix(iPoint, iDim, jDim));

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        for (size_t iDim = 0; iDim < nDim; ++iDim)
          AD::SetPreaccOut(gradient(iPoint, iVar, iDim));

      AD::EndPreacc();
    }

    /*--- Correct the gradient values across any periodic boundaries. ---*/

    if (solver != nullptr)
    {
      SU2_OMP_MASTER
      {
        for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
        {
          solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
          solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
        }
      }
      SU2_OMP_BARRIER
    }

    /*--- Second loop over points of the grid to compute final gradient. ---*/

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
    {
      /*--- Entries of upper triangular matrix R. ---*/

      su2double r11 = Rmatrix(iPoint,0,0);
      su2double r12 = Rmatrix(iPoint,0,1);
      su2double r22 = Rmatrix(iPoint,1,1);
      su2double r13 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;

      AD::StartPreacc();
      AD::SetPreaccIn(r11);
      AD::SetPreaccIn(r12);
      AD::SetPreaccIn(r22);

      if (r11 >= 0.0) r11 = sqrt(r11);
      if (r11 >= 0.0) r12 /= r11; else r12 = 0.0;
      su2double tmp = r22-r12*r12;
      if (tmp >= 0.0) r22 = sqrt(tmp); else r22 = 0.0;

      if (nDim == 3) {
        r13   = Rmatrix(iPoint,0,2);
        r23_a = Rmatrix(iPoint,1,2);
        r23_b = Rmatrix(iPoint,2,1);
        r33   = Rmatrix(iPoint,2,2);

        AD::SetPreaccIn(r13);
        AD::SetPreaccIn(r23_a);
        AD::SetPreaccIn(r23_b);
        AD::SetPreaccIn(r33);

        if (r11 >= 0.0) r13 /= r11; else r13 = 0.0;

        if ((r22 >= 0.0) && (r11*r22 >= 0.0)) {
          r23 = r23_a/r22 - r23_b*r12/(r11*r22);
        } else {
          r23 = 0.0;
        }

        tmp = r33 - r23*r23 - r13*r13;
        if (tmp >= 0.0) r33 = sqrt(tmp); else r33 = 0.0;
      }

      /*--- Compute determinant ---*/

      su2double detR2 = (r11*r22)*(r11*r22);
      if (nDim == 3) detR2 *= r33*r33;

      /*--- Detect singular matrices ---*/

      bool singular = false;

      if (detR2 <= EPS) {
        detR2 = 1.0;
        singular = true;
      }

      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/

      su2double Smatrix[MAXNDIM][MAXNDIM];

      if (singular) {
        for (size_t iDim = 0; iDim < nDim; ++iDim)
          for (size_t jDim = 0; jDim < nDim; ++jDim)
            Smatrix[iDim][jDim] = 0.0;
      }
      else {
        if (nDim == 2) {
          Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
          Smatrix[0][1] = -r11*r12/detR2;
          Smatrix[1][0] = Smatrix[0][1];
          Smatrix[1][1] = r11*r11/detR2;
        }
        else {
          su2double z11 = r22*r33;
          su2double z12 =-r12*r33;
          su2double z13 = r12*r23-r13*r22;
          su2double z22 = r11*r33;
          su2double z23 =-r11*r23;
          su2double z33 = r11*r22;

          Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
          Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
          Smatrix[0][2] = (z13*z33)/detR2;
          Smatrix[1][0] = Smatrix[0][1];
          Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
          Smatrix[1][2] = (z23*z33)/detR2;
          Smatrix[2][0] = Smatrix[0][2];
          Smatrix[2][1] = Smatrix[1][2];
          Smatrix[2][2] = (z33*z33)/detR2;
        }
      }

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        for (size_t jDim = 0; jDim < nDim; ++jDim)
          AD::SetPreaccOut(Smatrix[iDim][jDim]);

      AD::EndPreacc();

      /*--- Computation of the gradient: S*c ---*/

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        su2double Cvector[MAXNDIM];

        for (size_t iDim = 0; iDim < nDim; ++iDim)
        {
          Cvector[iDim] = 0.0;
          for (size_t jDim = 0; jDim < nDim; ++jDim)
            Cvector[iDim] += Smatrix[iDim][jDim] * gradient(iPoint, iVar, jDim);
        }

        for (size_t iDim = 0; iDim < nDim; ++iDim)
          gradient(iPoint, iVar, iDim) = Cvector[iDim];
      }
    }

  } // end SU2_OMP_PARALLEL

  /*--- If no solver was provided we do not communicate ---*/

  if (solver == nullptr) return;

  /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);

}
