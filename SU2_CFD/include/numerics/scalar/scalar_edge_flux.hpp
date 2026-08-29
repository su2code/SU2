/*!
 * \file scalar_edge_flux.hpp
 * \brief Model-agnostic convection and diffusion of a transported scalar, shared by every
 *        scalar solver (turbulence, transition, species, flamelet, heat).
 * \author P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/containers/container_decorators.hpp"
#include "../util.hpp"
#include "../../variables/CFlowVariable.hpp"

/*!
 * \brief Locates the data of one endpoint of an edge.
 * \note The j side of a boundary flux indexes per-marker ghost containers, which have the
 *       same types as the solver's own, so the kernels read both through one code path.
 */
template <class VariableType>
struct EdgeSide {
  const VariableType& scalarNodes;      /*!< \brief Scalar solver variables. */
  const CFlowVariable* flowNodes;       /*!< \brief Flow variables, null for solid heat transfer. */
  CMatrixView<const su2double> coord;   /*!< \brief Point coordinates. */
  CMatrixView<const su2double> gridVel; /*!< \brief Empty when the grid is static. */
};

/*!
 * \brief Loop invariant flags for a scalar edge flux, built once outside the edge loop so the
 *        compiler can unswitch the branches they guard.
 */
struct ScalarFluxOptions {
  bool dynamicGrid, boundedScalar, correctGradient, accurateJacobians;
  bool convective; /*!< \brief Whether the convective scheme contributes. */
  bool viscous;    /*!< \brief Whether the diffusion term contributes. */
  bool oneSided;   /*!< \brief Whether only the row of i is assembled. */
  bool muscl;      /*!< \brief Whether the convective scheme reconstructs. A boundary clears it. */
};

/*!
 * \brief Thin wrapper giving a Vector the all(iVar) accessor the MUSCL reconstruction helpers
 *        expect (see CCompressiblePrimitives in numerics_simd/flow/variables.hpp).
 */
template <class Double, size_t Size>
struct CScalarValues {
  static constexpr size_t nVar = Size; /*!< \brief Only used as VarType::nVar when reconstruct's own nVarGrad_ default
                                          applies; every call site here passes an explicit nVarGrad_ instead. */
  Vector<Double, Size> all;
};

/*!
 * \brief Diffusion of a transported scalar, driven by model-supplied coefficients.
 * \note The derived class returns the coefficients of both orientations of the edge, as one
 *       object, so that whatever the two share is computed once. A model whose coefficients
 *       are the same in both orientations returns the same value twice. A model whose matrix
 *       is diagonal declares DiagonalDiffusion and returns a vector instead of a matrix.
 */
template <class Double, class Derived, class FlowIndices, int nDim, size_t nVar>
class CAvgGradScalarBase {
 protected:
  using Int = typename CLaneTraits<Double>::Int;

  template <class VariableType>
  FORCEINLINE void diffusionTerms(const FlowIndices& idx, const ScalarFluxOptions& opt, Int iPoint,
                                  const EdgeSide<VariableType>& side_i, Int jPoint,
                                  const EdgeSide<VariableType>& side_j, const Vector<Double, nDim>& normal,
                                  const Vector<Double, nDim>& vector_ij, EdgeResidual<Double, nVar>& res) const {
    if (!opt.viscous) return;

    constexpr size_t Size = EdgeResidual<Double, nVar>::Size;

    const Double dist2_ij = fmax(squaredNorm(vector_ij), EPS);
    const Double proj_vector_ij = dot(vector_ij, normal) / dist2_ij;

    /*--- Average gradient, corrected for skewness when asked.
     * \note Gathered one variable at a time, bounded by res.nVar rather than Size: a static
     *       model with nVar 1 has Size 2 (the Matrix<Double,1,*> degeneracy floor), and a
     *       dynamic one has Size MaxScalarVar, so a single Size-wide read would run past the
     *       actual width of the gradient container in either case. ---*/
    Matrix<Double, Size, nDim> avgGrad;
    for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
      const auto grad_i = gatherVariables<1, nDim>(iPoint, side_i.scalarNodes.GetGradient(), iVar);
      const auto grad_j = gatherVariables<1, nDim>(jPoint, side_j.scalarNodes.GetGradient(), iVar);
      for (int iDim = 0; iDim < nDim; ++iDim) avgGrad(iVar, iDim) = 0.5 * (grad_i(iDim) + grad_j(iDim));
    }

    if (opt.correctGradient) {
      for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
        const Double phi_i = gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), iVar);
        const Double phi_j = gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), iVar);
        const Double corr = (dot(avgGrad[iVar], vector_ij) - phi_j + phi_i) / dist2_ij;
        for (int iDim = 0; iDim < nDim; ++iDim) avgGrad(iVar, iDim) -= corr * vector_ij(iDim);
      }
    }

    Vector<Double, Size> projGrad;
    for (size_t iVar = 0; iVar < res.nVar; ++iVar) projGrad(iVar) = dot(avgGrad[iVar], normal);

    /*--- The Jacobians of a conservative model are w.r.t. the conserved (density-weighted)
     * variable, which divides the geometric projection by the density of the row being written. ---*/
    Double w_i = 1.0, w_j = 1.0;
    if constexpr (Derived::Conservative) {
      w_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
      w_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());
    }
    const Double proj_on_w_i = proj_vector_ij / w_i;
    const Double proj_on_w_j = proj_vector_ij / w_j;

    const auto* self = static_cast<const Derived*>(this);
    const auto D = self->coefficients(idx, iPoint, side_i, jPoint, side_j);

    for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
      if constexpr (Derived::DiagonalDiffusion) {
        res.flux_i(iVar) -= D.i(iVar) * projGrad(iVar);
        res.jac_ii(iVar, iVar) += D.i(iVar) * proj_on_w_i;
        res.jac_ij(iVar, iVar) -= D.i(iVar) * proj_on_w_j;

        if (!opt.oneSided) {
          res.flux_j(iVar) += D.j(iVar) * projGrad(iVar);
          res.jac_ji(iVar, iVar) -= D.j(iVar) * proj_on_w_i;
          res.jac_jj(iVar, iVar) += D.j(iVar) * proj_on_w_j;
        }
      } else {
        for (size_t jVar = 0; jVar < res.nVar; ++jVar) {
          res.flux_i(iVar) -= D.i(iVar, jVar) * projGrad(jVar);
          res.jac_ii(iVar, jVar) += D.i(iVar, jVar) * proj_on_w_i;
          res.jac_ij(iVar, jVar) -= D.i(iVar, jVar) * proj_on_w_j;

          if (!opt.oneSided) {
            res.flux_j(iVar) += D.j(iVar, jVar) * projGrad(jVar);
            res.jac_ji(iVar, jVar) -= D.j(iVar, jVar) * proj_on_w_i;
            res.jac_jj(iVar, jVar) += D.j(iVar, jVar) * proj_on_w_j;
          }
        }
      }
    }

    if (opt.accurateJacobians) {
      /*--- Coefficients that depend on the transported variables contribute here. A model whose
       * correction is a per-edge constant (e.g. SA's) can ignore the side/point arguments; one
       * whose correction depends on point values (e.g. SST's, on the transported variable at
       * either endpoint) needs them, so every model is handed the same full context diffusionTerms
       * itself has, matching extraDiffusionTerms's signature below. ---*/
      self->coefficientJacobians(idx, iPoint, side_i, jPoint, side_j, projGrad, res);
    }

    self->extraDiffusionTerms(idx, iPoint, side_i, jPoint, side_j, normal, vector_ij, res);
  }

  /*!
   * \brief Contribution of the derivatives of the coefficients themselves.
   */
  template <class... Ts>
  FORCEINLINE void coefficientJacobians(Ts&...) const {}

  /*!
   * \brief Diffusion of a model that transports more than one gradient, of states it
   *        synthesises from its own containers.
   */
  template <class... Ts>
  FORCEINLINE void extraDiffusionTerms(Ts&...) const {}
};

/*!
 * \brief Convective flux shared by every model whose transport equation has the shape
 *        flux(iVar) = a0 * w_i * phi_i(iVar) + a1 * w_j * phi_j(iVar), which is every
 *        model except SA and stochastic backscatter (see CScalarFlux_SA).
 * \note The weight w is 1 for a non-conservative model, the density for a conservative one.
 */
template <class Double, class Derived, class FlowIndices, int nDim, size_t nVar>
class CUpwScalarFlux : public CAvgGradScalarBase<Double, Derived, FlowIndices, nDim, nVar> {
 protected:
  using Int = typename CLaneTraits<Double>::Int;

  explicit CUpwScalarFlux(const CConfig&) {}

  /*!
   * \param[in] phi - Transported variable of both endpoints, reconstructed if opt.muscl is set;
   *            read from here rather than side_i/side_j.scalarNodes directly so a model needs
   *            no reconstruction logic of its own.
   */
  template <class VariableType, size_t Size>
  FORCEINLINE void finalizeFlux(const FlowIndices& idx, const ScalarFluxOptions&, Int iPoint,
                                const EdgeSide<VariableType>& side_i, Int jPoint, const EdgeSide<VariableType>& side_j,
                                const Double& a0, const Double& a1, const CPair<CScalarValues<Double, Size>>& phi,
                                EdgeResidual<Double, nVar>& res) const {
    Double w0 = a0, w1 = a1;
    if constexpr (Derived::Conservative) {
      w0 *= gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
      w1 *= gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());
    }

    for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
      const Double flux = w0 * phi.i.all(iVar) + w1 * phi.j.all(iVar);

      res.flux_i(iVar) += flux;
      res.flux_j(iVar) -= flux;

      res.jac_ii(iVar, iVar) += w0;
      res.jac_ij(iVar, iVar) += w1;
      res.jac_ji(iVar, iVar) -= w0;
      res.jac_jj(iVar, iVar) -= w1;
    }
  }
};

/*!
 * \brief Upwind convection and diffusion of a transported scalar, accumulated into one
 *        residual, each term contributing or not according to the options.
 */
template <class Double_, class Derived, class FlowIndices, int nDim_, size_t nVar_>
class CUpwScalarBase : public CUpwScalarFlux<Double_, Derived, FlowIndices, nDim_, nVar_> {
 public:
  using Double = Double_;
  using Int = typename CLaneTraits<Double_>::Int;
  static constexpr int nDim = nDim_;
  static constexpr size_t nVar = nVar_;

  /*!
   * \brief Backing size, for the model to size the containers its coefficients() returns;
   *        nVar itself is Dynamic for a runtime model and never usable as a container size.
   */
  static constexpr size_t Size = EdgeResidual<Double_, nVar_>::Size;

 protected:
  using Base = CUpwScalarFlux<Double_, Derived, FlowIndices, nDim_, nVar_>;

  const FlowIndices idx;
  const size_t nEqn; /*!< \brief Equations of the model, which a dynamic one gives to its base. */

  /*!
   * \brief MUSCL reconstruction parameters, read from CConfig once per construction (i.e. once
   *        per nonlinear iteration, see CScalarSolver::EdgeFluxResidual) instead of per edge;
   *        this is also where the scalar limiter's freezing (GetLimiterIter) is resolved, by
   *        collapsing its type to NONE once frozen. The flow limiter is not frozen this way: once
   *        the flow solver stops recomputing it, it keeps applying the last values it has.
   */
  const su2double kappa, umusclRamp, kappaFlow;
  const LIMITER limiterType, limiterTypeFlow;
  const bool musclFlow;

 public:
  /*!
   * \brief Constructor, inherited by the model with `using Base::Base`.
   * \note Public, not protected: a using-declaration that inherits a constructor keeps the
   *       base's own access, so the solver that builds the concrete model needs this public
   *       to build it at all.
   */
  explicit CUpwScalarBase(const CConfig& config, size_t nEqn_ = nVar_)
      : Base(config),
        idx(nDim, config.GetnSpecies()),
        nEqn(nEqn_),
        kappa(config.GetMUSCL_Kappa()),
        umusclRamp(config.GetMUSCLRampValue()),
        kappaFlow(config.GetMUSCL_Kappa_Flow()),
        limiterType(config.GetInnerIter() <= config.GetLimiterIter() ? config.GetKind_SlopeLimit() : LIMITER::NONE),
        limiterTypeFlow(config.GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE ? config.GetKind_SlopeLimit_Flow()
                                                                                     : LIMITER::NONE),
        musclFlow(config.GetMUSCL_Flow() && config.GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) {}

  template <class VariableType>
  FORCEINLINE EdgeResidual<Double, nVar> ComputeFlux(const ScalarFluxOptions& opt, Int iPoint,
                                                     const EdgeSide<VariableType>& side_i, Int jPoint,
                                                     const EdgeSide<VariableType>& side_j,
                                                     const Vector<Double, nDim>& normal, const Double& massFlux) const {
    /*--- Inputs are registered as they are read, by each of the two terms. ---*/
    AD::StartPreacc();
    AD::SetPreaccIn(normal, nDim);

    EdgeResidual<Double, nVar> res(nEqn);

    /*--- Read once by the reconstruction and by the diffusion. ---*/
    Vector<Double, nDim> vector_ij;
    if (opt.muscl || opt.viscous) {
      vector_ij = distanceVector<nDim>(iPoint, side_i.coord, jPoint, side_j.coord);
    }

    if (opt.convective) {
      /*--- Upwinding weights of the face normal mass or volume flux. ---*/
      Double a0, a1;
      if (opt.boundedScalar) {
        AD::SetPreaccIn(massFlux);
        const Double rho_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
        const Double rho_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());
        a0 = fmax(0.0, massFlux) / rho_i;
        a1 = fmin(0.0, massFlux) / rho_j;
      } else {
        /*--- The mass-flux branch above reads the edge flux computed from unreconstructed flow
         * primitives directly, so only this branch needs a reconstructed velocity. ---*/
        CPair<CScalarValues<Double, nDim>> u;
        u.i.all = gatherVariables<nDim>(iPoint, side_i.flowNodes->GetPrimitive(), idx.Velocity());
        u.j.all = gatherVariables<nDim>(jPoint, side_j.flowNodes->GetPrimitive(), idx.Velocity());

        if (opt.muscl && musclFlow) {
          reconstruct<nDim>(iPoint, jPoint, vector_ij, side_i.flowNodes->GetGradient_Reconstruction(),
                            side_i.flowNodes->GetLimiter_Primitive(), limiterTypeFlow, idx.Velocity(), u, kappaFlow,
                            umusclRamp);
        }

        /*--- Face normal velocity of the mean of the two points, relative to the grid. ---*/
        Vector<Double, nDim> vel_ij;
        for (int iDim = 0; iDim < nDim; ++iDim) vel_ij(iDim) = 0.5 * (u.i.all(iDim) + u.j.all(iDim));

        if (opt.dynamicGrid) {
          const auto ug_i = gatherVariables<nDim>(iPoint, side_i.gridVel);
          /*--- A boundary's ghost point has no grid velocity of its own: it is spatially
           * coincident with i, so it moves with it. ---*/
          const auto ug_j = opt.oneSided ? ug_i : gatherVariables<nDim>(jPoint, side_j.gridVel);
          for (int iDim = 0; iDim < nDim; ++iDim) vel_ij(iDim) -= 0.5 * (ug_i(iDim) + ug_j(iDim));
        }

        const Double q_ij = dot(vel_ij, normal);
        a0 = fmax(0.0, q_ij);
        a1 = fmin(0.0, q_ij);
      }

      /*--- Transported variable of both endpoints, reconstructed if opt.muscl is set.
       * Gathered one variable at a time (like the diffusion gradients above) so a static model
       * with nVar 1 never reads past the single column its solution container actually has. ---*/
      CPair<CScalarValues<Double, Size>> phi;
      for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
        phi.i.all(iVar) = gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), iVar);
        phi.j.all(iVar) = gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), iVar);
      }

      if constexpr (nVar != Dynamic) {
        if (opt.muscl) {
          reconstruct<nVar>(iPoint, jPoint, vector_ij, side_i.scalarNodes.GetGradient_Reconstruction(),
                            side_i.scalarNodes.GetLimiter(), limiterType, 0, phi, kappa, umusclRamp);
        }
      }

      static_cast<const Derived*>(this)->finalizeFlux(idx, opt, iPoint, side_i, jPoint, side_j, a0, a1, phi, res);
    }

    Base::diffusionTerms(idx, opt, iPoint, side_i, jPoint, side_j, normal, vector_ij, res);

    AD::SetPreaccOut(res.flux_i, res.nVar);
    if (!opt.oneSided) AD::SetPreaccOut(res.flux_j, res.nVar);
    AD::EndPreacc();

    return res;
  }

  /*!
   * \brief Compute the flux of an edge and write it to the linear system.
   */
  template <class VariableType>
  FORCEINLINE void ComputeFlux(const ScalarFluxOptions& opt, Int iEdge, Int iPoint,
                               const EdgeSide<VariableType>& side_i, Int jPoint, const EdgeSide<VariableType>& side_j,
                               const Vector<Double, nDim>& normal, const Double& massFlux, bool implicit,
                               UpdateType updateType, Double updateMask, CSysVector<su2double>& vector,
                               CSysVector<su2double>& vectorDiff, SparseMatrixType& matrix) const {
    const auto res = ComputeFlux(opt, iPoint, side_i, jPoint, side_j, normal, massFlux);

    updateLinearSystem(iEdge, iPoint, jPoint, implicit, updateType, updateMask, res, vector, vectorDiff, matrix);
  }
};
