/*!
 * \file computeGradientsGreenGauss.hpp
 * \brief Generic implementation of Green-Gauss gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
 * \author P. Gomes
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include <vector>
#include <algorithm>

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "computeGradientsSymmetry.hpp"

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
template <size_t nDim, class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver, MPI_QUANTITIES kindMpiComm, PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry, const CConfig& config, const FieldType& field, size_t varBegin,
                                size_t varEnd, GradientType& gradient) {
  const size_t nPointDomain = geometry.GetnPointDomain();
  bool isFlowSolver = solver->GetSolverName().find("FLOW") != string::npos;

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  static constexpr size_t MAXNVAR = 20;
  static constexpr size_t MAXNDIM = 3;
  static constexpr size_t MAXNSYMS = 100;

  /*--- For each (non-halo) volume integrate over its faces (edges). ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    auto nodes = geometry.nodes;

    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();
    AD::SetPreaccIn(nodes->GetVolume(iPoint));
    AD::SetPreaccIn(nodes->GetPeriodicVolume(iPoint));

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar) AD::SetPreaccIn(field(iPoint, iVar));

    /*--- Clear the gradient. --*/

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim) gradient(iPoint, iVar, iDim) = 0.0;

    /*--- Handle averaging and division by volume in one constant. ---*/

    su2double halfOnVol = 0.5 / (nodes->GetVolume(iPoint) + nodes->GetPeriodicVolume(iPoint));

    /*--- Add a contribution due to each neighbor. ---*/

    for (size_t iNeigh = 0; iNeigh < nodes->GetnPoint(iPoint); ++iNeigh) {
      size_t iEdge = nodes->GetEdge(iPoint, iNeigh);
      size_t jPoint = nodes->GetPoint(iPoint, iNeigh);

      /*--- Determine if edge points inwards or outwards of iPoint.
       *    If inwards we need to flip the area vector. ---*/

      su2double dir = (iPoint < jPoint) ? 1.0 : -1.0;
      su2double weight = dir * halfOnVol;

      const auto area = geometry.edges->GetNormal(iEdge);
      AD::SetPreaccIn(area, nDim);

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
        AD::SetPreaccIn(field(jPoint, iVar));
        su2double flux = weight * (field(iPoint, iVar) + field(jPoint, iVar));

        for (size_t iDim = 0; iDim < nDim; ++iDim) gradient(iPoint, iVar, iDim) += flux * area[iDim];
      }
    }

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim) AD::SetPreaccOut(gradient(iPoint, iVar, iDim));

    AD::EndPreacc();
  }
  END_SU2_OMP_FOR

  su2double flux[MAXNVAR] = {0.0};

  /*--- Add edges of markers that contribute to the gradients ---*/
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if ((config.GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config.GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
        (config.GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE) &&
        (config.GetMarker_All_KindBC(iMarker) != EULER_WALL) &&
        (config.GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

      /*--- Work is shared in inner loop as two markers
       *    may try to update the same point. ---*/

      SU2_OMP_FOR_STAT(32)
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
        auto nodes = geometry.nodes;

        /*--- Halo points do not need to be considered. ---*/

        if (!nodes->GetDomain(iPoint)) continue;

        su2double volume = nodes->GetVolume(iPoint) + nodes->GetPeriodicVolume(iPoint);
        const auto area = geometry.vertex[iMarker][iVertex]->GetNormal();

        for (size_t iVar = varBegin; iVar < varEnd; iVar++)
          flux[iVar] = field(iPoint,iVar) / volume;

        for (size_t iVar = varBegin; iVar < varEnd; iVar++) {
          for (size_t iDim = 0; iDim < nDim; iDim++) {
            gradient(iPoint, iVar, iDim) -= flux[iVar] * area[iDim];
          }
        } // loop over variables
      } // vertices
      END_SU2_OMP_FOR
    } //found right marker
  } // iMarkers

  /* For symmetry planes, we need to impose the conditions (Blazek eq. 8.40):
   * 1. n.grad(phi) = 0
   * 2. n.grad(v.t) = 0
   * 3. t.grad(v.n) = 0
   */

  /*--- Check how many symmetry planes there are ---*/
  unsigned short Syms[MAXNSYMS] = {0};
  unsigned short nSym = 0;
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if ((config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) ||
       (config.GetMarker_All_KindBC(iMarker) == EULER_WALL)) {
    Syms[nSym] = iMarker;
    nSym++;
    }
  }

  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {

    if ((config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) || (config.GetMarker_All_KindBC(iMarker) == EULER_WALL)) {
      SU2_OMP_FOR_STAT(32)
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {

        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();

        /*--- Normal vector for this vertex (negate for outward convention). ---*/
        const su2double* VertexNormal = geometry.vertex[iMarker][iVertex]->GetNormal();

        const auto NormArea = GeometryToolbox::Norm(nDim, VertexNormal);

        su2double UnitNormal[nDim] = {0.0};
        for (size_t iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = VertexNormal[iDim] / NormArea;

        /*--- Normal of the primary symmetry plane ---*/
        su2double NormalPrim[MAXNDIM] = {0.0}, UnitNormalPrim[MAXNDIM] = {0.0};

        /*--- At this point we can find out if the node is shared with another symmetry.
         * Step 1: do we have other symmetries? ---*/
        if (nSym>1) {
          /*--- Step 2: are we on a shared node? ---*/
          for (auto jMarker=0;jMarker<nSym;jMarker++) {
            /*--- we do not need the current symmetry ---*/
            if (iMarker!= Syms[jMarker]) {
              /*--- Loop over all points on the other symmetry and check if current iPoint is on the symmetry ---*/
              for (auto jVertex = 0ul; jVertex < geometry.nVertex[Syms[jMarker]]; jVertex++) {
                const auto jPoint = geometry.vertex[Syms[jMarker]][jVertex]->GetNode();
                if (iPoint==jPoint) {
                  /*--- Does the other symmetry have a lower ID? Then that is the primary symmetry ---*/
                  if (Syms[jMarker]<iMarker) {
                    /*--- So whe have to get the normal of that other marker ---*/
                    geometry.vertex[Syms[jMarker]][jVertex]->GetNormal(NormalPrim);
                    su2double AreaPrim = GeometryToolbox::Norm(nDim, NormalPrim);
                    for(unsigned short iDim = 0; iDim < nDim; iDim++) {
                      UnitNormalPrim[iDim] = NormalPrim[iDim] / AreaPrim;
                    }

                    /*--- Correct the current normal as n2_new = n2 - (n2.n1)n1 ---*/
                    su2double ProjNorm = 0.0;
                    for (auto iDim = 0u; iDim < nDim; iDim++) ProjNorm += UnitNormal[iDim] * UnitNormalPrim[iDim];
                    /*--- We check if the normal of the 2 planes coincide.
                     * We only update the normal if the normals of the symmetry planes are different. ---*/
                    if (fabs(1.0-ProjNorm)>EPS) {
                      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] -= ProjNorm * UnitNormalPrim[iDim];
                      /* Make normalized vector ---*/
                      su2double newarea=GeometryToolbox::Norm(nDim, UnitNormal);
                      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = UnitNormal[iDim]/newarea;
                    }
                  }
                }
              }
            }
          }
        }


        /*--- Tensor mapping from global Cartesian to local Cartesian aligned with symmetry plane ---*/
        su2activematrix TensorMap(nDim,nDim);

        /*--- Compute a new base for TensorMap aligned with the unit normal ---*/
        BaseFromNormal(nDim,UnitNormal,TensorMap);

        su2activematrix Gradients_iPoint(varEnd-varBegin,nDim);

        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            Gradients_iPoint[iVar][iDim] = gradient(iPoint, iVar, iDim);
          }
        }

        ReflectGradient(nDim, varBegin,varEnd, isFlowSolver, TensorMap, Gradients_iPoint);

        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            gradient(iPoint, iVar, iDim) = Gradients_iPoint[iVar][iDim];
          }
        }

      } // loop over vertices
      END_SU2_OMP_FOR
    } // symmetry
  } // markers

  /*--- If no solver was provided we do not communicate ---*/

  if (solver == nullptr) return;

  /*--- Account for periodic contributions. ---*/

  for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic() / 2; ++iPeriodic) {
    solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
  }

  /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);
}
}  // namespace detail



/*!
 * \brief Instantiations for 2D and 3D.
 * \ingroup FvmAlgos
 */
template <class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver, MPI_QUANTITIES kindMpiComm, PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry, const CConfig& config, const FieldType& field, size_t varBegin,
                                size_t varEnd, GradientType& gradient) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::computeGradientsGreenGauss<2>(solver, kindMpiComm, kindPeriodicComm, geometry, config, field, varBegin,
                                            varEnd, gradient);
      break;
    case 3:
      detail::computeGradientsGreenGauss<3>(solver, kindMpiComm, kindPeriodicComm, geometry, config, field, varBegin,
                                            varEnd, gradient);
      break;
    default:
      SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
      break;
  }
}
