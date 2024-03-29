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
 bool isFlowSolver = solver->GetSolverName().find("FLOW") != string::npos;

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  static constexpr size_t MAXNVAR = 20;
  static constexpr size_t MAXNDIM = 3;
  static constexpr size_t MAXNSYMS = 5;
  //bool sym2 = false;

  /*--- Allocation of primitive gradient arrays for viscous fluxes. ---*/
  su2activematrix Grad_Reflected(varEnd, nDim);

  /*--- Tensor mapping from global Cartesian to local Cartesian aligned with symmetry plane ---*/
  su2activematrix TensorMap(nDim,nDim);
  su2activematrix TensorMapInv(nDim,nDim);
  su2activematrix Gradients_Velocity(nDim,nDim);
  su2activematrix Gradients_Velocity_Reflected(nDim,nDim);

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

  /*--- Check how many symmetry planes there are ---*/
  unsigned short Syms[MAXNSYMS] = {0};
  unsigned short nSym = 0;
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) {
    Syms[nSym] = iMarker;
    nSym++;
    }
  }

  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) {
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {

        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();

        /*--- Normal vector for this vertex (negate for outward convention). ---*/
        const su2double* VertexNormal = geometry.vertex[iMarker][iVertex]->GetNormal();

        const auto NormArea = GeometryToolbox::Norm(nDim, VertexNormal);

        su2double UnitNormal[nDim] = {0.0};
        for (size_t iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = VertexNormal[iDim] / NormArea;

  // normal of the primary symmetry plane
    su2double NormalPrim[MAXNDIM] = {0.0}, UnitNormalPrim[MAXNDIM] = {0.0};

    // at this point we can find out if the node is shared with another symmetry.
    // step 1: do we have other symmetries?
    if (nSym>1) {
      // step 2: are we on a shared node?
      for (auto jMarker=0;jMarker<nSym;jMarker++) {
        // we do not need the current symmetry
        if (iMarker!= Syms[jMarker]) {
          // loop over all points on the other symmetry and check if current iPoint is on the symmetry

          for (auto jVertex = 0ul; jVertex < geometry.nVertex[Syms[jMarker]]; jVertex++) {
            const auto jPoint = geometry.vertex[Syms[jMarker]][jVertex]->GetNode();

            const su2double *coor = geometry.nodes->GetCoord(jPoint);

            if (iPoint==jPoint) {
              // Does the other symmetry have a lower ID? Then that is the primary symmetry
              if (Syms[jMarker]<iMarker) {
                // so whe have to get the normal of that other marker
                geometry.vertex[Syms[jMarker]][jVertex]->GetNormal(NormalPrim);
                su2double AreaPrim = GeometryToolbox::Norm(nDim, NormalPrim);
                for(unsigned short iDim = 0; iDim < nDim; iDim++) {
                  UnitNormalPrim[iDim] = NormalPrim[iDim] / AreaPrim;
                }

               // correct the current normal as n2_new = n2 - (n2.n1).n1
                su2double ProjNorm = 0.0;
                for (auto iDim = 0u; iDim < nDim; iDim++) ProjNorm += UnitNormal[iDim] * UnitNormalPrim[iDim];
                // We check if the normal of the 2 planes coincide. We only update the normal if the normals of the symmetry planes are different.
                if (fabs(1.0-ProjNorm)>EPS) {
                  for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] -= ProjNorm * UnitNormalPrim[iDim];
                  // make normalized vector again
                  su2double newarea=GeometryToolbox::Norm(nDim, UnitNormal);
                  for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = UnitNormal[iDim]/newarea;
                }

              }
            }

          }
        }
      }
    }

    if (isFlowSolver == true) {

        /*--- Preprocessing: Compute unit tangential, the direction is arbitrary as long as
              t*n=0 && |t|_2 = 1 ---*/
        su2double TangentialNorm, Tangential[MAXNDIM] = {0.0};
        su2double Orthogonal[MAXNDIM] = {0.0};
        switch (nDim) {
          case 2: {
            Tangential[0] = -UnitNormal[1];
            Tangential[1] = UnitNormal[0];
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              TensorMap[0][iDim] = UnitNormal[iDim];
              TensorMap[1][iDim] = Tangential[iDim];
              //TensorMap[2][iDim] = Orthogonal[iDim];
              TensorMapInv[iDim][0] = UnitNormal[iDim];
              TensorMapInv[iDim][1] = Tangential[iDim];
              //TensorMapInv[iDim][2] = Orthogonal[iDim];
            }
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

            /*--- Compute 3rd direction of the base using cross product ---*/
            Orthogonal[0] = UnitNormal[1]*Tangential[2] - UnitNormal[2]*Tangential[1];
            Orthogonal[1] = UnitNormal[2]*Tangential[0] - UnitNormal[0]*Tangential[2];
            Orthogonal[2] = UnitNormal[0]*Tangential[1] - UnitNormal[1]*Tangential[0];

            // now we construct the tensor mapping T, note that its inverse is the transpose of T
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              TensorMap[0][iDim] = UnitNormal[iDim];
              TensorMap[1][iDim] = Tangential[iDim];
              TensorMap[2][iDim] = Orthogonal[iDim];
              TensorMapInv[iDim][0] = UnitNormal[iDim];
              TensorMapInv[iDim][1] = Tangential[iDim];
              TensorMapInv[iDim][2] = Orthogonal[iDim];
            }
            break;
          }
        }  // switch


        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            Grad_Reflected[iVar][iDim] = gradient(iPoint, iVar, iDim);
          }
        }


        /*--- Get gradients of primitives of boundary cell ---*/
        for (auto iVar = 0; iVar < nDim; iVar++)
          for (auto iDim = 0u; iDim < nDim; iDim++)
            Gradients_Velocity[iVar][iDim] = gradient(iPoint, iVar+1, iDim);


        for (auto iDim = 0u; iDim < nDim; iDim++) {
          for (auto jDim = 0u; jDim < nDim; jDim++) {
            Gradients_Velocity_Reflected[iDim][jDim] = 0.0;
          }
        }

        // Q' = L^T*Q*T
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          for (auto jDim = 0u; jDim < nDim; jDim++) {
            for (auto kDim = 0u; kDim < nDim; kDim++) {
              for (auto mDim = 0u; mDim < nDim; mDim++) {
                Gradients_Velocity_Reflected[iDim][jDim] += TensorMap[iDim][mDim]*TensorMap[jDim][kDim]*Gradients_Velocity[mDim][kDim];
              }
            }
          }
        }


        // now also remove the gradients that are supposed to be zero
        //first entry is normal direction n, with velocity V,W in the sym plane
        // and velocity U in the normal direction
        // we need to remove dV/dn and DW/dn = 0
        // and dU/dt and dU/ds = 0

        // we have aligned such that U is the direction of the normal
        // in 2D: dU/dy = dV/dx = 0
        // in 3D: dU/dy = dV/dx = 0
        //        dU/dz = dW/dx = 0
        for (auto iDim = 1u; iDim < nDim; iDim++) {
          Gradients_Velocity_Reflected[0][iDim] = 0.0;
          Gradients_Velocity_Reflected[iDim][0] = 0.0;
        }

       for (auto iDim = 0u; iDim < nDim; iDim++) {
          for (auto jDim = 0u; jDim < nDim; jDim++) {
            Gradients_Velocity[iDim][jDim] = 0.0;
          }
        }

        // now transform back by taking the inverse again
       // T = (L^-1)*T'
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          for (auto jDim = 0u; jDim < nDim; jDim++) {
            for (auto kDim = 0u; kDim < nDim; kDim++) {
              for (auto mDim = 0u; mDim < nDim; mDim++) {
                Gradients_Velocity[iDim][jDim] += TensorMap[mDim][iDim]*TensorMap[kDim][jDim]*Gradients_Velocity_Reflected[mDim][kDim];
              }
            }
          }
        }

       for (auto iDim = 0u; iDim < nDim; iDim++) {
           for (auto jDim = 0u; jDim < nDim; jDim++) {
             Grad_Reflected[iDim+1][jDim] = Gradients_Velocity[iDim][jDim];
           }
         }
    }

        /*--- Reflect the gradients for all scalars including the velocity components.
              The gradients of the velocity components are set later with the
              correct values: grad(V)_r = grad(V) - 2 [grad(V)*n]n, V being any primitive ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          // For the auxiliary variable we do not have velocity vectors. But the gradients of the auxvars
          // do not seem to be used so this has no effect on the computations.
          if (iVar == 0 || iVar > nDim || isFlowSolver==true) {  // We should exclude here for instance AuxVars for axisymmetry?

            /*--- Compute projected part of the gradient in a dot product ---*/
            su2double ProjGradient = 0.0;
            for (auto iDim = 0u; iDim < nDim; iDim++) ProjGradient += Grad_Reflected[iVar][iDim] * UnitNormal[iDim];

            for (auto iDim = 0u; iDim < nDim; iDim++)
              // full reflection
              //Grad_Reflected[iVar][iDim] = Grad_Reflected[iVar][iDim] - 2.0 * ProjGradient * UnitNormal[iDim];
              // reflect once to eliminate normal direction
              Grad_Reflected[iVar][iDim] = Grad_Reflected[iVar][iDim] - ProjGradient * UnitNormal[iDim];
          }
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
