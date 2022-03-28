/*!
 * \file computeGradientsL2Projection.hpp
 * \brief Generic implementation of L2 projection gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
 * \author B. Munguía
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

namespace detail {

/*!
 * \brief Compute the gradient of a field using the L2 projection theorem.
 * \note Template nDim to allow efficient unrolling of inner loops.
 * \note Gradients can be computed only for a contiguous range of variables, defined
 *       by [varBegin, varEnd[ (e.g. 0,1 computes the gradient of the 1st variable).
 *       This can be used, for example, to compute only velocity gradients.
 * \note The function uses an optional solver object to perform communications, if
 *       none (nullptr) is provided the function does not fail (the objective of
 *       this is to improve test-ability).
 * \note Only implemented for tri (2D) and tet (3D) elements.
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
void computeGradientsL2Projection(CSolver* solver,
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
  const size_t nElem = geometry.GetnElem();
  const size_t nFace = (nDim == 2) ? 3 : 4;
  const size_t nNode = nFace - 1; // Number of nodes on a face

  su2double* normal = new su2double[nDim];

  const su2double factor = (nDim == 2) ? 1.0/6.0 : 1.0/24.0;

  /*--- For each (non-halo) volume integrate over its faces (edges). ---*/

  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
  {
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = 0.0;
  }
  
  /*--- For each volume integrate over its faces (edges). ---*/
  auto nodes = geometry.nodes;
  for (size_t iElem = 0; iElem < nElem; ++iElem)
  {
    auto elem = geometry.elem[iElem];

    /*--- Add the contribution from each node of the volume. ---*/
    for (size_t iFace = 0; iFace < nFace; ++iFace) {
      const size_t iNode = nNode - iFace;
      const size_t iPoint = elem->GetNode(iNode);

      /*--- Inward normal of opposite face ---*/
      if (nDim == 2) {
        const size_t jPoint = elem->GetNode(elem->GetFaces(iFace,0));
        const size_t kPoint = elem->GetNode(elem->GetFaces(iFace,1));
        normal[0] = nodes->GetCoord(jPoint, 1) - nodes->GetCoord(kPoint, 1);
        normal[1] = nodes->GetCoord(kPoint, 0) - nodes->GetCoord(jPoint, 0);
      }
      /*---TODO: Implement cross-product for tets ---*/
      else {
      }

      /*--- Gradient contribution of the node ---*/
      for (size_t jNode = 0; jNode < nNode; ++jNode) {
        const size_t jPoint = elem->GetNode(jNode);
        if (!nodes->GetDomain(jPoint)) continue;
        const su2double Vol = nodes->GetVolume(jPoint);
        for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
          const su2double var = field(iPoint,iVar);
          for (size_t iDim = 0; iDim < nDim; ++iDim) {
            gradient(jPoint, iVar, iDim) += factor*var*normal[iDim]/Vol;
          }
        }
      }
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

template<size_t nDim, class GradientType>
void computeHessiansL2Projection(CSolver* solver,
                                 MPI_QUANTITIES kindMpiComm,
                                 PERIODIC_QUANTITIES kindPeriodicComm,
                                 CGeometry& geometry,
                                 const CConfig& config,
                                 const GradientType& gradient,
                                 size_t varBegin,
                                 size_t varEnd,
                                 GradientType& hessian)
{
  const size_t nPointDomain = geometry.GetnPointDomain();
  const size_t nElem = geometry.GetnElem();
  const size_t nFace = (nDim == 2) ? 3 : 4;
  const size_t nNode = nFace - 1; // Number of nodes on a face

  su2double* normal = new su2double[nDim];

  const su2double factor = (nDim == 2) ? 1.0/6.0 : 1.0/24.0;

  /*--- For each (non-halo) volume integrate over its faces (edges). ---*/

  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
  {
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < 3*(nDim-1); ++iDim)
        hessian(iPoint, iVar, iDim) = 0.0;
  }
  
  /*--- For each volume integrate over its faces (edges). ---*/
  auto nodes = geometry.nodes;
  for (size_t iElem = 0; iElem < nElem; ++iElem)
  {
    auto elem = geometry.elem[iElem];

    /*--- Add the contribution from each node of the volume. ---*/
    for (size_t iFace = 0; iFace < nFace; ++iFace) {
      const size_t iNode = nNode - iFace;
      const size_t iPoint = elem->GetNode(iNode);

      /*--- Inward normal of opposite face ---*/
      if (nDim == 2) {
        const size_t jPoint = elem->GetNode(elem->GetFaces(iFace,0));
        const size_t kPoint = elem->GetNode(elem->GetFaces(iFace,1));
        normal[0] = nodes->GetCoord(jPoint, 1) - nodes->GetCoord(kPoint, 1);
        normal[1] = nodes->GetCoord(kPoint, 0) - nodes->GetCoord(jPoint, 0);
      }
      /*---TODO: Implement cross-product for tets ---*/
      else {
      }

      /*--- Gradient contribution of the node ---*/
      for (size_t jNode = 0; jNode < nNode; ++jNode) {
        const size_t jPoint = elem->GetNode(jNode);
        if (!nodes->GetDomain(jPoint)) continue;
        const su2double Vol = nodes->GetVolume(jPoint);
        for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
          for (size_t iDim = 0; iDim < nDim; ++iDim) {
            for (size_t jDim = 0; jDim < nDim; ++jDim) {
              su2double grad = gradient(iPoint,iVar,jDim);
              size_t ind = (iDim <= jDim) ? iDim*nDim - ((iDim - 1)*iDim)/2 + jDim - iDim 
                                          : jDim*nDim - ((jDim - 1)*jDim)/2 + iDim - jDim;
              if (iDim != jDim) grad *= 0.5;
              hessian(jPoint, iVar, ind) += factor*grad*normal[iDim]/Vol;
            }
          }
        }
      }
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
 */
template<class FieldType, class GradientType>
void computeGradientsL2Projection(CSolver* solver,
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
    detail::computeGradientsL2Projection<2>(solver, kindMpiComm, kindPeriodicComm, geometry,
                                            config, field, varBegin, varEnd, gradient);
    break;
  case 3:
    detail::computeGradientsL2Projection<3>(solver, kindMpiComm, kindPeriodicComm, geometry,
                                            config, field, varBegin, varEnd, gradient);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}

template<class GradientType>
void computeHessiansL2Projection(CSolver* solver,
                                 MPI_QUANTITIES kindMpiComm,
                                 PERIODIC_QUANTITIES kindPeriodicComm,
                                 CGeometry& geometry,
                                 const CConfig& config,
                                 const GradientType& gradient,
                                 size_t varBegin,
                                 size_t varEnd,
                                 GradientType& hessian) {
  switch (geometry.GetnDim()) {
  case 2:
    detail::computeHessiansL2Projection<2>(solver, kindMpiComm, kindPeriodicComm, geometry,
                                           config, gradient, varBegin, varEnd, hessian);
    break;
  case 3:
    detail::computeHessiansL2Projection<3>(solver, kindMpiComm, kindPeriodicComm, geometry,
                                           config, gradient, varBegin, varEnd, hessian);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}
