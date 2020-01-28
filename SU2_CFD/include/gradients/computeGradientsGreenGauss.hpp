/*!
 * \file computeGradientsGreenGauss.hpp
 * \brief Generic implementation of Green-Gauss gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
 * \author P. Gomes
 * \version 7.0.0 "Blackbird"
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
 * \brief Compute the gradient of a field using the Green-Gauss theorem.
 * \note Gradients can be computed only for a contiguous range of variables, defined
 *       by [varBegin, varEnd[ (e.g. 0,1 computes the gradient of the 1st variable).
 *       This can be used, for example, to compute only velocity gradients.
 * \note The function uses an optional solver object to perform communications, if
 *       none (nullptr) is provided the function does not fail (the objective of
 *       this is to improve test-ability).
 * \param[in] solver - Optional, solver associated with the field (used only for MPI).
 * \param[in] kindMpiComm - Type of MPI communication required.
 * \param[in] kindPeriodicComm - Type of periodic communication required.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem, used to identify types of boundaries.
 * \param[in] field - Generic object implementing operator (iPoint, iVar).
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 */
template<class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver,
                                MPI_QUANTITIES kindMpiComm,
                                PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry,
                                CConfig& config,
                                const FieldType& field,
                                size_t varBegin,
                                size_t varEnd,
                                GradientType& gradient)
{
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
    /*--- For each (non-halo) volume integrate over its faces (edges). ---*/

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
    {
      auto node = geometry.node[iPoint];

      AD::StartPreacc();
      AD::SetPreaccIn(node->GetVolume());
      AD::SetPreaccIn(node->GetPeriodicVolume());

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        AD::SetPreaccIn(field(iPoint,iVar));

      /*--- Clear the gradient. --*/

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        for (size_t iDim = 0; iDim < nDim; ++iDim)
          gradient(iPoint, iVar, iDim) = 0.0;

      /*--- Handle averaging and division by volume in one constant. ---*/

      su2double halfOnVol = 0.5 / (node->GetVolume()+node->GetPeriodicVolume());

      /*--- Add a contribution due to each neighbor. ---*/

      for (size_t iNeigh = 0; iNeigh < node->GetnPoint(); ++iNeigh)
      {
        size_t iEdge = node->GetEdge(iNeigh);
        size_t jPoint = node->GetPoint(iNeigh);

        /*--- Determine if edge points inwards or outwards of iPoint.
         *    If inwards we need to flip the area vector. ---*/

        su2double dir = (iPoint == geometry.edge[iEdge]->GetNode(0))? 1.0 : -1.0;
        su2double weight = dir * halfOnVol;

        const su2double* area = geometry.edge[iEdge]->GetNormal();
        AD::SetPreaccIn(area, nDim);

        for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        {
          AD::SetPreaccIn(field(jPoint,iVar));

          su2double flux = weight * (field(iPoint,iVar) + field(jPoint,iVar));

          for (size_t iDim = 0; iDim < nDim; ++iDim)
            gradient(iPoint, iVar, iDim) += flux * area[iDim];
        }

      }

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        for (size_t iDim = 0; iDim < nDim; ++iDim)
          AD::SetPreaccOut(gradient(iPoint,iVar,iDim));

      AD::EndPreacc();
    }

    /*--- Add boundary fluxes. ---*/

    for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker)
    {
      if ((config.GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config.GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY))
      {
        /*--- Work is shared in inner loop as two markers
         *    may try to update the same point. ---*/

        SU2_OMP_FOR_STAT(OMP_MAX_CHUNK)
        for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex)
        {
          size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
          auto node = geometry.node[iPoint];

          /*--- Halo points do not need to be considered. ---*/

          if (!node->GetDomain()) continue;

          su2double volume = node->GetVolume() + node->GetPeriodicVolume();

          const su2double* area = geometry.vertex[iMarker][iVertex]->GetNormal();

          for (size_t iVar = varBegin; iVar < varEnd; iVar++)
          {
            su2double flux = field(iPoint,iVar) / volume;

            for (size_t iDim = 0; iDim < nDim; iDim++)
              gradient(iPoint, iVar, iDim) -= flux * area[iDim];
          }
        }
      }
    }

  } // end SU2_OMP_PARALLEL

  /*--- If no solver was provided we do not communicate ---*/

  if (solver == nullptr) return;

  /*--- Account for periodic contributions. ---*/

  for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
  {
    solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
  }

  /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);

}
