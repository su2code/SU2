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
 * \note Built through the named constructors below: the flags are too many and too alike to be
 *       given positionally, where one transposed pair would change the discretization silently.
 */
struct ScalarFluxOptions {
  bool dynamicGrid = false;
  bool boundedScalar = false;
  bool correctGradient = false;
  bool accurateJacobians = false;
  bool implicit = false;  /*!< \brief Whether the Jacobians are assembled, and so computed. */
  bool convective = true; /*!< \brief Whether the convective scheme contributes. */
  bool viscous = false;   /*!< \brief Whether the diffusion term contributes. */
  bool oneSided = false;  /*!< \brief Whether only the row of i is assembled. */
  bool muscl = false;     /*!< \brief Whether the convective scheme reconstructs. */

  /*!
   * \brief Options of the interior edge loop: both terms, both rows, and reconstruction and
   *        gradient correction as the configuration asks for them.
   */
  static ScalarFluxOptions Interior(const CConfig& config, bool bounded, bool accurateJacobians = false) {
    auto opt = common(config);
    opt.boundedScalar = bounded;
    opt.accurateJacobians = accurateJacobians;
    opt.correctGradient = true;
    opt.viscous = true;
    opt.muscl = config.GetMUSCL();
    return opt;
  }

  /*!
   * \brief Options of a boundary that imposes a convective flux alone, which is most of them:
   *        the diffusive term at an inlet or an outlet causes serious convergence problems.
   */
  static ScalarFluxOptions BoundaryConvective(const CConfig& config, bool bounded) {
    auto opt = common(config);
    opt.boundedScalar = bounded;
    opt.oneSided = true;
    return opt;
  }

  /*!
   * \brief Options of a boundary that imposes both terms, which is the turbomachinery sites.
   * \note They impose no mass flux, so the bounded scheme contributes nothing here whatever the
   *       configuration says.
   */
  static ScalarFluxOptions BoundaryFull(const CConfig& config) {
    auto opt = common(config);
    opt.correctGradient = true;
    opt.viscous = true;
    opt.oneSided = true;
    return opt;
  }

  /*!
   * \brief Options of the diffusive pass of a fluid interface, which follows a convective pass
   *        over the donor vertices.
   * \param[in] correctGrad - Whether the projected gradient is corrected for skewness, which the
   *            models do not agree on at this boundary.
   */
  static ScalarFluxOptions BoundaryDiffusive(const CConfig& config, bool correctGrad) {
    auto opt = common(config);
    opt.correctGradient = correctGrad;
    opt.convective = false;
    opt.viscous = true;
    opt.oneSided = true;
    return opt;
  }

 private:
  static ScalarFluxOptions common(const CConfig& config) {
    ScalarFluxOptions opt;
    opt.dynamicGrid = config.GetDynamic_Grid();
    opt.implicit = config.GetKind_TimeIntScheme() == EULER_IMPLICIT;
    return opt;
  }
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

  /*!
   * \param[in] rho - Density of both endpoints, read once by ComputeFlux.
   */
  template <class VariableType>
  FORCEINLINE void diffusionTerms(const FlowIndices& idx, const ScalarFluxOptions& opt, Int iPoint,
                                  const EdgeSide<VariableType>& side_i, Int jPoint,
                                  const EdgeSide<VariableType>& side_j, const CPair<Double>& rho,
                                  const Vector<Double, nDim>& normal, const Vector<Double, nDim>& vector_ij,
                                  EdgeResidual<Double, nVar>& res) const {
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
      const auto grad_i = gatherVariables<nDim>(iPoint, side_i.scalarNodes.GetGradient(), iVar);
      const auto grad_j = gatherVariables<nDim>(jPoint, side_j.scalarNodes.GetGradient(), iVar);
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
    Double proj_on_w_i = proj_vector_ij, proj_on_w_j = proj_vector_ij;
    if constexpr (Derived::Conservative) {
      proj_on_w_i = proj_vector_ij / rho.i;
      proj_on_w_j = proj_vector_ij / rho.j;
    }

    const auto* self = static_cast<const Derived*>(this);
    const auto D = self->coefficients(idx, iPoint, side_i, jPoint, side_j, rho);

    for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
      if constexpr (Derived::DiagonalDiffusion) {
        res.flux_i(iVar) -= D.i(iVar) * projGrad(iVar);
        if (opt.implicit) {
          res.jac_ii(iVar, iVar) += D.i(iVar) * proj_on_w_i;
          if (!opt.oneSided) res.jac_ij(iVar, iVar) -= D.i(iVar) * proj_on_w_j;
        }
        if (!opt.oneSided) {
          res.flux_j(iVar) += D.j(iVar) * projGrad(iVar);
          if (opt.implicit) {
            res.jac_ji(iVar, iVar) -= D.j(iVar) * proj_on_w_i;
            res.jac_jj(iVar, iVar) += D.j(iVar) * proj_on_w_j;
          }
        }
      } else {
        for (size_t jVar = 0; jVar < res.nVar; ++jVar) {
          res.flux_i(iVar) -= D.i(iVar, jVar) * projGrad(jVar);
          if (opt.implicit) {
            res.jac_ii(iVar, jVar) += D.i(iVar, jVar) * proj_on_w_i;
            if (!opt.oneSided) res.jac_ij(iVar, jVar) -= D.i(iVar, jVar) * proj_on_w_j;
          }
          if (!opt.oneSided) {
            res.flux_j(iVar) += D.j(iVar, jVar) * projGrad(jVar);
            if (opt.implicit) {
              res.jac_ji(iVar, jVar) -= D.j(iVar, jVar) * proj_on_w_i;
              res.jac_jj(iVar, jVar) += D.j(iVar, jVar) * proj_on_w_j;
            }
          }
        }
      }
    }

    if (opt.implicit && opt.accurateJacobians) {
      /*--- Coefficients that depend on the transported variables contribute here, from whatever
       * the model chose to carry in the object it returned from coefficients. ---*/
      self->coefficientJacobians(opt, D, projGrad, res);
    }

    self->extraDiffusionTerms(idx, opt, iPoint, side_i, jPoint, side_j, rho, normal, vector_ij, res);
  }

  /*!
   * \brief Contribution of the derivatives of the coefficients themselves.
   */
  template <class... Ts>
  FORCEINLINE void coefficientJacobians(Ts&&...) const {}

  /*!
   * \brief Diffusion of a model that transports more than one gradient, of states it
   *        synthesises from its own containers.
   */
  template <class... Ts>
  FORCEINLINE void extraDiffusionTerms(Ts&&...) const {}
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
   * \brief Upwind convection of the transported variable, weighted by the density for a
   *        conservative model.
   * \param[in] phi - Transported variable of both endpoints, reconstructed if opt.muscl is set;
   *            read from here rather than side_i/side_j.scalarNodes directly so a model needs
   *            no reconstruction logic of its own.
   * \param[in] rho - Density of both endpoints, reconstructed alongside the velocity when the
   *            convective scheme reconstructs, so it weights the flux as the velocity does.
   * \note The flux is written in terms of the transported variable but the Jacobians are w.r.t.
   *       the conserved one, which for a conservative model is the density-weighted variable;
   *       the density therefore multiplies the flux and not the Jacobian.
   */
  template <class VariableType, size_t Size>
  FORCEINLINE void finalizeFlux(const FlowIndices&, const ScalarFluxOptions& opt, Int, const EdgeSide<VariableType>&,
                                Int, const EdgeSide<VariableType>&, const Double& a0, const Double& a1,
                                const CPair<Double>& rho, const CPair<CScalarValues<Double, Size>>& phi,
                                EdgeResidual<Double, nVar>& res) const {
    Double w0 = a0, w1 = a1;
    if constexpr (Derived::Conservative) {
      w0 *= rho.i;
      w1 *= rho.j;
    }

    for (size_t iVar = 0; iVar < res.nVar; ++iVar) {
      const Double flux = w0 * phi.i.all(iVar) + w1 * phi.j.all(iVar);

      res.flux_i(iVar) += flux;
      if (!opt.oneSided) res.flux_j(iVar) -= flux;

      if (opt.implicit) {
        res.jac_ii(iVar, iVar) += a0;
        if (!opt.oneSided) {
          res.jac_ij(iVar, iVar) += a1;
          res.jac_ji(iVar, iVar) -= a0;
          res.jac_jj(iVar, iVar) -= a1;
        }
      }
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

  /*!
   * \brief Whether the model's diffusion coefficients read the density, beyond the reading that
   *        Conservative already implies. A model that declares neither never gathers it.
   */
  static constexpr bool DiffusionReadsDensity = false;

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
        musclFlow(config.GetMUSCL_Flow() && config.GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) {
    if (nEqn > Size) {
      SU2_MPI::Error("Static arrays are too small for the requested equation count.", CURRENT_FUNCTION);
    }
  }

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

    /*--- Density of both endpoints, gathered once: the conservative weighting of the convective
     * term, the bounded scheme's division of the mass flux, and some models' diffusion
     * coefficients all want it, and in reverse mode every gather is a preaccumulation input. ---*/
    CPair<Double> rho{Double(1.0), Double(1.0)};
    if (Derived::Conservative || Derived::DiffusionReadsDensity || opt.boundedScalar) {
      rho.i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
      rho.j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());
    }

    if (opt.convective) {
      /*--- Upwinding weights of the face normal mass or volume flux, and the density that weights
       * a conservative flux, which follows the velocity in being reconstructed or not. ---*/
      Double a0, a1;
      CPair<Double> rhoConv = rho;

      if (opt.boundedScalar) {
        AD::SetPreaccIn(massFlux);
        a0 = fmax(0.0, massFlux) / rho.i;
        a1 = fmin(0.0, massFlux) / rho.j;
      } else {
        CPair<CScalarValues<Double, nDim>> u;
        u.i.all = gatherVariables<nDim>(iPoint, side_i.flowNodes->GetPrimitive(), idx.Velocity());
        u.j.all = gatherVariables<nDim>(jPoint, side_j.flowNodes->GetPrimitive(), idx.Velocity());

        if (opt.muscl && musclFlow) {
          reconstruct<nDim>(iPoint, jPoint, vector_ij, side_i.flowNodes->GetGradient_Reconstruction(),
                            side_i.flowNodes->GetLimiter_Primitive(), limiterTypeFlow, idx.Velocity(), u, kappaFlow,
                            umusclRamp);

          if constexpr (Derived::Conservative) {
            /*--- Density is not adjacent to the velocity in the primitive row, so it is a second
             * reconstruction of one variable rather than a wider slice of the first. ---*/
            CPair<CScalarValues<Double, 1>> r;
            r.i.all(0) = rho.i;
            r.j.all(0) = rho.j;
            reconstruct<1>(iPoint, jPoint, vector_ij, side_i.flowNodes->GetGradient_Reconstruction(),
                           side_i.flowNodes->GetLimiter_Primitive(), limiterTypeFlow, idx.Density(), r, kappaFlow,
                           umusclRamp);
            rhoConv.i = r.i.all(0);
            rhoConv.j = r.j.all(0);
          }
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

      if (opt.muscl) {
        if constexpr (nVar != Dynamic) {
          reconstruct<nVar>(iPoint, jPoint, vector_ij, side_i.scalarNodes.GetGradient_Reconstruction(),
                            side_i.scalarNodes.GetLimiter(), limiterType, 0, phi, kappa, umusclRamp);
        } else {
          /*--- A dynamic model's equation count is only known at runtime, so the reconstructed
           * width is passed as an argument instead of a template parameter. ---*/
          reconstruct(iPoint, jPoint, vector_ij, side_i.scalarNodes.GetGradient_Reconstruction(),
                      side_i.scalarNodes.GetLimiter(), limiterType, 0, phi, kappa, umusclRamp, res.nVar);
        }
      }

      static_cast<const Derived*>(this)->finalizeFlux(idx, opt, iPoint, side_i, jPoint, side_j, a0, a1, rhoConv, phi,
                                                      res);
    }

    Base::diffusionTerms(idx, opt, iPoint, side_i, jPoint, side_j, rho, normal, vector_ij, res);

    setPreaccOut(res.flux_i, res.nVar);
    if (!opt.oneSided) setPreaccOut(res.flux_j, res.nVar);
    AD::EndPreacc();

    return res;
  }

  /*!
   * \brief Compute the flux of an edge and write it to the linear system.
   */
  template <class VariableType>
  FORCEINLINE void ComputeFlux(const ScalarFluxOptions& opt, Int iEdge, Int iPoint,
                               const EdgeSide<VariableType>& side_i, Int jPoint, const EdgeSide<VariableType>& side_j,
                               const Vector<Double, nDim>& normal, const Double& massFlux, UpdateType updateType,
                               Double updateMask, CSysVector<su2double>& vector, CSysVector<su2double>& vectorDiff,
                               SparseMatrixType& matrix) const {
    const auto res = ComputeFlux(opt, iPoint, side_i, jPoint, side_j, normal, massFlux);

    updateLinearSystem(iEdge, iPoint, jPoint, opt.implicit, updateType, updateMask, res, vector, vectorDiff, matrix);
  }
};
