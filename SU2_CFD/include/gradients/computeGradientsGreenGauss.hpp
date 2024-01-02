/*!
 * \file computeGradientsGreenGauss.hpp
 * \brief Generic implementation of Green-Gauss gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
 * \author P. Gomes
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
 * \brief Compute the gradient of a field using the Green-Gauss theorem.
 * \ingroup FvmAlgos
 * \note Template nDim to allow efficient unrolling of inner loops.
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
template<size_t nDim, class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver,
                                MPI_QUANTITIES kindMpiComm,
                                PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry,
                                const CConfig& config,
                                const FieldType& field,
                                size_t varBegin,
                                size_t varEnd,
                                GradientType& gradient)
{
  const size_t nPointDomain = geometry.GetnPointDomain();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- For each (non-halo) volume integrate over its faces (edges). ---*/


  // ********************************************************************
  // loop over all cells
  // ********************************************************************
  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
  {
    auto nodes = geometry.nodes;

    if (iPoint==255){
      cout << "iPoint = "<<iPoint<<endl;
    }

    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();
    AD::SetPreaccIn(nodes->GetVolume(iPoint));
    AD::SetPreaccIn(nodes->GetPeriodicVolume(iPoint));

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      AD::SetPreaccIn(field(iPoint,iVar));

    /*--- Clear the gradient. --*/

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = 0.0;

    /*--- Handle averaging and division by volume in one constant. ---*/

    su2double halfOnVol = 0.5 / (nodes->GetVolume(iPoint)+nodes->GetPeriodicVolume(iPoint));

    /*--- Add a contribution due to each neighbor. ---*/

    for (size_t iNeigh = 0; iNeigh < nodes->GetnPoint(iPoint); ++iNeigh)
    {
      size_t iEdge = nodes->GetEdge(iPoint,iNeigh);
      size_t jPoint = nodes->GetPoint(iPoint,iNeigh);

      /*--- Determine if edge points inwards or outwards of iPoint.
       *    If inwards we need to flip the area vector. ---*/

      su2double dir = (iPoint < jPoint)? 1.0 : -1.0;
      su2double weight = dir * halfOnVol;

      const auto area = geometry.edges->GetNormal(iEdge);
      AD::SetPreaccIn(area, nDim);

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        AD::SetPreaccIn(field(jPoint,iVar));
        su2double vali = field(iPoint,iVar);
        su2double valj = field(jPoint,iVar);
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
  END_SU2_OMP_FOR


  // ********************************************************************
  // loop over all cells on a symmetry plane
  // ********************************************************************

  /*--- Add GG boundary fluxes on symmetry. ---*/
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker)
  {
    if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) {
      //cout << "symmetry plane found" <<endl;
      SU2_OMP_FOR_STAT(32)
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
        // recompute gradients on symmetry plane
        auto nodes = geometry.nodes;
        /*--- Clear the gradient. --*/
        for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
          for (size_t iDim = 0; iDim < nDim; ++iDim)
            gradient(iPoint, iVar, iDim) = 0.0;
        /*--- Handle averaging and division by volume in one constant. ---*/
        su2double halfOnVol = 0.5 / (nodes->GetVolume(iPoint)+nodes->GetPeriodicVolume(iPoint));
        /*--- Add a contribution due to each neighbor. ---*/
        for (size_t iNeigh = 0; iNeigh < nodes->GetnPoint(iPoint); ++iNeigh)
        {
          size_t iEdge = nodes->GetEdge(iPoint,iNeigh);
          size_t jPoint = nodes->GetPoint(iPoint,iNeigh);
          su2double *coordi={nodes->GetCoord(iPoint)};
          su2double *coordj={nodes->GetCoord(jPoint)};
          /*--- Determine if edge points inwards or outwards of iPoint.
           *    If inwards we need to flip the area vector. ---*/
          su2double dir = (iPoint < jPoint)? 1.0 : -1.0;
          su2double weight = dir * halfOnVol;
          const auto area = geometry.edges->GetNormal(iEdge);
          for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
          {
            su2double vali = field(iPoint,iVar);
            su2double valj = field(jPoint,iVar);
            su2double flux = weight * (field(iPoint,iVar) + field(jPoint,iVar));
            for (size_t iDim = 0; iDim < nDim; ++iDim)
              gradient(iPoint, iVar, iDim) += flux * area[iDim];
          }
        }
      }
      END_SU2_OMP_FOR
    }
  }

  // ********************************************************************
  // loop over all cells on other boundaries
  // ********************************************************************
  /*--- Add boundary fluxes. ---*/
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker)
  {
    if ((config.GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config.GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
        (config.GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE) &&
        (config.GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY))
    {

      /*--- Work is shared in inner loop as two markers
       *    may try to update the same point. ---*/

      SU2_OMP_FOR_STAT(32)
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex)
      {
        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
        auto nodes = geometry.nodes;

        /*--- Halo points do not need to be considered. ---*/

        if (!nodes->GetDomain(iPoint)) continue;

        su2double volume = nodes->GetVolume(iPoint) + nodes->GetPeriodicVolume(iPoint);

        const su2double* area = geometry.vertex[iMarker][iVertex]->GetNormal();

        for (size_t iVar = varBegin; iVar < varEnd; iVar++)
        {
          su2double flux = field(iPoint,iVar) / volume;

          for (size_t iDim = 0; iDim < nDim; iDim++)
            gradient(iPoint, iVar, iDim) -= flux * area[iDim];
        }
      }
      END_SU2_OMP_FOR
    }
  }



  /*--- Add boundary flux correction for symmetry. ---*/
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker)
  {
      if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) {
        //cout << "symmetry plane found" <<endl;
        SU2_OMP_FOR_STAT(32)
        for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
          size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
          // recompute gradients on symmetry plane
          auto nodes = geometry.nodes;
          /*--- Halo points do not need to be considered. ---*/
          if (!nodes->GetDomain(iPoint)) continue;
          /*--- get the reflected state [phi_x,phi_y,phi_z]---*/
          su2double Area;
          su2double value;
          su2double gradphi[3] = {0.0};
          su2double phi[3] = {0.0};
          su2double phi_reflected[3]={0.0};
          su2double Normal[3]={0.0};
          su2double UnitNormal[3]={0.0};
          // now we construct field_x,field_y from field(iPoint,iVar)
          /*--- Normal vector for this vertex (negate for outward convention). ---*/
          geometry.vertex[iMarker][iVertex]->GetNormal(Normal);
          for (size_t iDim = 0; iDim < nDim; iDim++) 
            Normal[iDim] = -Normal[iDim];
          Area = GeometryToolbox::Norm(nDim, Normal);
          for (size_t iDim = 0; iDim < nDim; iDim++) 
            UnitNormal[iDim] = -Normal[iDim] / Area;

          // loop over all iVar and set gradient normal to the boundary to zero
          // But for velocity, the velocity gradient normal to the wall is not zero,
          // only the velocity gradient tangential to the wall is zero
          for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
            value = field(iPoint,iVar);
            for (size_t iDim = 0; iDim < nDim; iDim++)
              gradphi[iDim] = gradient(iPoint, iVar, iDim);
            // we can now project the vector gradphi to the wall-aligned coordinates
            // compute grad(phi).n
            su2double Projphi = 0.0;
            for (size_t iDim = 0; iDim < nDim; iDim++)
              Projphi += gradphi[iDim]*UnitNormal[iDim];
            // mirror reflection
            for (size_t iDim = 0; iDim < nDim; iDim++)
              // complete reflection:
              //phi_reflected[iDim] = gradphi[iDim] - 2.0 * Projphi * UnitNormal[iDim];
              // tangential component only:
              // this can only work when the wall is aligned with the cartesian coordinates
              //if (iVar==2){
              //  // for the velocity normal to the symmetry, the tangential gradient is zero.
              //  phi_reflected[iDim] = Projphi * UnitNormal[iDim];
             //  
             // }
             // else{
                // all scalars and the tangential velocity gradients have normal component zero (only tangential component left)
                phi_reflected[iDim] = gradphi[iDim] - 1.0 * Projphi * UnitNormal[iDim];
             // }
  
            //cout << "modify gradient" << endl;
            for (size_t iDim = 0; iDim < nDim; iDim++)
              gradient(iPoint, iVar, iDim) = phi_reflected[iDim];
          }
        }
      END_SU2_OMP_FOR
      }
  }



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
} // end namespace

/*!
 * \brief Instantiations for 2D and 3D.
 * \ingroup FvmAlgos
 */
template<class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver,
                                MPI_QUANTITIES kindMpiComm,
                                PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry,
                                const CConfig& config,
                                const FieldType& field,
                                size_t varBegin,
                                size_t varEnd,
                                GradientType& gradient) {
  switch (geometry.GetnDim()) {
  case 2:
    detail::computeGradientsGreenGauss<2>(solver, kindMpiComm, kindPeriodicComm, geometry,
                                          config, field, varBegin, varEnd, gradient);
    break;
  case 3:
    detail::computeGradientsGreenGauss<3>(solver, kindMpiComm, kindPeriodicComm, geometry,
                                          config, field, varBegin, varEnd, gradient);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}
