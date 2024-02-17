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

#include <vector>
#include <algorithm>

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

namespace detail {

// find local vertex on a symmetry marker using global iPoint
inline su2double* getVertexNormalfromPoint(const CConfig& config, CGeometry& geometry, unsigned long iPointGlobal){
  unsigned long iPointSym=0;
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) {
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
        iPointSym = geometry.vertex[iMarker][iVertex]->GetNode();
        if (iPointSym == iPointGlobal)
          return geometry.vertex[iMarker][iVertex]->GetNormal();
      }
    }
  }
  cout << "point is not found " << endl;
  exit(0);
}

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


cout << "Green Gauss: solver name = " << solver->GetSolverName() << endl;
cout << "number of variables = " << varEnd << endl;
cout << "commtype= = " << kindMpiComm << endl;
cout << "viscous = " << config.GetViscous();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  static constexpr size_t MAXNVAR = 20;
  static constexpr size_t MAXNDIM = 3;

  /*--- Allocation of primitive gradient arrays for viscous fluxes. ---*/
  su2activematrix Grad_Reflected(varEnd, nDim);

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

  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) {
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {

        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
        auto nodes = geometry.nodes;
        // we need to set the gradient to zero for the entire marker to prevent double-counting
        // points that are shared by other markers
        //for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        //  for (size_t iDim = 0; iDim < nDim; ++iDim) gradient(iPoint, iVar, iDim) = 0.0;

        su2double halfOnVol = 0.5 / (nodes->GetVolume(iPoint) + nodes->GetPeriodicVolume(iPoint));

        /*--- Normal vector for this vertex (negate for outward convention). ---*/
        const su2double* VertexNormal = geometry.vertex[iMarker][iVertex]->GetNormal();

        // reflected normal V = U - 2*U_t
        const auto NormArea = GeometryToolbox::Norm(nDim, VertexNormal);

        su2double UnitNormal[nDim] = {0.0};
        for (size_t iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = VertexNormal[iDim] / NormArea;

        /*--- Preprocessing: Compute unit tangential, the direction is arbitrary as long as
              t*n=0 && |t|_2 = 1 ---*/
        su2double TangentialNorm, Tangential[MAXNDIM] = {0.0};
        switch (nDim) {
          case 2: {
            Tangential[0] = -UnitNormal[1];
            Tangential[1] = UnitNormal[0];
            break;
          }
          case 3: {
            /*--- n = ai + bj + ck, if |b| > |c| ---*/
            if (abs(UnitNormal[1]) > abs(UnitNormal[2])) {
              /*--- t = bi + (c-a)j - bk  ---*/
              Tangential[0] = UnitNormal[1];
              Tangential[1] = UnitNormal[2] - UnitNormal[0];
              Tangential[2] = -UnitNormal[1];
            } else {
              /*--- t = ci - cj + (b-a)k  ---*/
              Tangential[0] = UnitNormal[2];
              Tangential[1] = -UnitNormal[2];
              Tangential[2] = UnitNormal[1] - UnitNormal[0];
            }
            /*--- Make it a unit vector. ---*/
            TangentialNorm = sqrt(pow(Tangential[0], 2) + pow(Tangential[1], 2) + pow(Tangential[2], 2));
            Tangential[0] = Tangential[0] / TangentialNorm;
            Tangential[1] = Tangential[1] / TangentialNorm;
            Tangential[2] = Tangential[2] / TangentialNorm;
            break;
          }
        }  // switch


        /*--- Get gradients of primitives of boundary cell ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++)
          for (auto iDim = 0u; iDim < nDim; iDim++)
            Grad_Reflected[iVar][iDim] = gradient(iPoint, iVar, iDim);

        /*--- Reflect the gradients for all scalars including the velocity components.
              The gradients of the velocity components are set later with the
              correct values: grad(V)_r = grad(V) - 2 [grad(V)*n]n, V being any primitive ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          // For the auxiliary variable we do not have velocity vectors. But the gradients of the auxvars
          // do not seem to be used so this has no effect on the computations.
          if (iVar == 0 || iVar > nDim) {  // We should exclude here for instance AuxVars for axisymmetry?

            /*--- Compute projected part of the gradient in a dot product ---*/
            su2double ProjGradient = 0.0;
            for (auto iDim = 0u; iDim < nDim; iDim++) ProjGradient += Grad_Reflected[iVar][iDim] * UnitNormal[iDim];

            for (auto iDim = 0u; iDim < nDim; iDim++)
              Grad_Reflected[iVar][iDim] = Grad_Reflected[iVar][iDim] - 2.0 * ProjGradient * UnitNormal[iDim];
          }
        }

        /*--- Compute gradients of normal and tangential velocity:
              grad(v*n) = grad(v_x) n_x + grad(v_y) n_y (+ grad(v_z) n_z)
              grad(v*t) = grad(v_x) t_x + grad(v_y) t_y (+ grad(v_z) t_z) ---*/

        /*--- if we do not have auxiliary gradients ---*/
        su2double GradNormVel[MAXNVAR] = {0.0};
        su2double GradTangVel[MAXNVAR] = {0.0};
        for (auto iVar = 0u; iVar < nDim; iVar++) {  // counts gradient components
          GradNormVel[iVar] = 0.0;
          GradTangVel[iVar] = 0.0;
          for (auto iDim = 0u; iDim < nDim; iDim++) {  // counts sum with unit normal/tangential
            GradNormVel[iVar] += Grad_Reflected[iDim + 1][iVar] * UnitNormal[iDim];
            GradTangVel[iVar] += Grad_Reflected[iDim + 1][iVar] * Tangential[iDim];
          }
        }

        /*--- Reflect gradients in tangential and normal direction by substracting the normal/tangential
              component twice, just as done with velocity above.
              grad(v*n)_r = grad(v*n) - 2 {grad([v*n])*t}t
              grad(v*t)_r = grad(v*t) - 2 {grad([v*t])*n}n ---*/
        su2double ProjNormVelGrad = 0.0;
        su2double ProjTangVelGrad = 0.0;
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          ProjNormVelGrad += GradNormVel[iDim] * Tangential[iDim];  // grad([v*n])*t
          ProjTangVelGrad += GradTangVel[iDim] * UnitNormal[iDim];  // grad([v*t])*n
        }

        for (auto iDim = 0u; iDim < nDim; iDim++) {
          GradNormVel[iDim] = GradNormVel[iDim] - 2.0 * ProjNormVelGrad * Tangential[iDim];
          GradTangVel[iDim] = GradTangVel[iDim] - 2.0 * ProjTangVelGrad * UnitNormal[iDim];
        }

        /*--- Transfer reflected gradients back into the Cartesian Coordinate system:
              grad(v_x)_r = grad(v*n)_r n_x + grad(v*t)_r t_x
              grad(v_y)_r = grad(v*n)_r n_y + grad(v*t)_r t_y
              ( grad(v_z)_r = grad(v*n)_r n_z + grad(v*t)_r t_z ) ---*/
        for (auto iVar = 0u; iVar < nDim; iVar++)    // loops over the velocity component gradients
          for (auto iDim = 0u; iDim < nDim; iDim++) {  // loops over the entries of the above
            Grad_Reflected[iVar + 1][iDim] =
                GradNormVel[iDim] * UnitNormal[iVar] + GradTangVel[iDim] * Tangential[iVar];
                // at this point, the gradients are done.
                //cout << "grad = " << Grad_Reflected[iVar + 1][iDim] << " " << gradient(iPoint,iVar+1,iDim) << endl;
            }

        /*--- Update gradients with reflected gradients ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++)
          for (auto iDim = 0u; iDim < nDim; iDim++)
            gradient(iPoint,iVar,iDim) = Grad_Reflected[iVar][iDim];

      } //ivertex
    } //symmetry
  } //loop over markers


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
