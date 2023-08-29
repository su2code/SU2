/*!
 * \file scalar_convection.hpp
 * \brief Declarations of numerics classes for discretization of
 *        convective fluxes in scalar problems.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "../CNumerics.hpp"

/*!
 * \class CUpwScalar
 * \brief Template class for scalar upwind fluxes between nodes i and j.
 * \details This class serves as a template for the scalar upwinding residual
 *   classes.  The general structure of a scalar upwinding calculation is the
 *   same for many different  models, which leads to a lot of repeated code.
 *   By using the template design pattern, these sections of repeated code are
 *   moved to this shared base class, and the specifics of each model
 *   are implemented by derived classes.  In order to add a new residual
 *   calculation for a convection residual, extend this class and implement
 *   the pure virtual functions with model-specific behavior.
 * \ingroup ConvDiscr
 * \author C. Pederson, A. Bueno., and A. Campos.
 */
template <class FlowIndices>
class CUpwScalar : public CNumerics {
 protected:
  enum : unsigned short {MAXNVAR = 8};

  const FlowIndices idx;            /*!< \brief Object to manage the access to the flow primitives. */
  su2double a0 = 0.0;               /*!< \brief The maximum of the face-normal velocity and 0. */
  su2double a1 = 0.0;               /*!< \brief The minimum of the face-normal velocity and 0. */
  su2double Flux[MAXNVAR];          /*!< \brief Final result, diffusive flux/residual. */
  su2double* Jacobian_i[MAXNVAR];   /*!< \brief Flux Jacobian w.r.t. node i. */
  su2double* Jacobian_j[MAXNVAR];   /*!< \brief Flux Jacobian w.r.t. node j. */
  su2double JacobianBuffer[2*MAXNVAR*MAXNVAR];  /*!< \brief Static storage for the two Jacobians. */

  const bool incompressible = false, dynamic_grid = false;

  /*!
   * \brief A pure virtual function. Derived classes must use it to register the additional
   *        variables they use as preaccumulation inputs, e.g. the density for SST.
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method, derived classes
   *        compute the Flux and its Jacobians via this method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(const CConfig* config) = 0;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwScalar(unsigned short ndim, unsigned short nvar, const CConfig* config)
    : CNumerics(ndim, nvar, config),
      idx(ndim, config->GetnSpecies()),
      incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE),
      dynamic_grid(config->GetDynamic_Grid()) {
    if (nVar > MAXNVAR) {
      SU2_MPI::Error("Static arrays are too small.", CURRENT_FUNCTION);
    }
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = &JacobianBuffer[iVar * nVar];
      Jacobian_j[iVar] = &JacobianBuffer[iVar * nVar + MAXNVAR * MAXNVAR];
    }

    /*--- Initialize the JacobianBuffer to zero. ---*/
    for (unsigned short iVar = 0; iVar < 2*MAXNVAR*MAXNVAR; iVar++) {
      JacobianBuffer[iVar] = 0.0;
    }
  }

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  CNumerics::ResidualType<> ComputeResidual(const CConfig* config) final {
    AD::StartPreacc();
    AD::SetPreaccIn(Normal, nDim);
    AD::SetPreaccIn(ScalarVar_i, nVar);
    AD::SetPreaccIn(ScalarVar_j, nVar);
    if (dynamic_grid) {
      AD::SetPreaccIn(GridVel_i, nDim);
      AD::SetPreaccIn(GridVel_j, nDim);
    }
    AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
    AD::SetPreaccIn(&V_j[idx.Velocity()], nDim);
    AD::SetPreaccIn(V_i[idx.Density()]);
    AD::SetPreaccIn(V_j[idx.Density()]);
    AD::SetPreaccIn(MassFlux);

    ExtraADPreaccIn();

    if (bounded_scalar) {
      a0 = fmax(0.0, MassFlux) / V_i[idx.Density()];
      a1 = fmin(0.0, MassFlux) / V_j[idx.Density()];
    } else {
      su2double q_ij = 0.0;
      if (dynamic_grid) {
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          su2double Velocity_i = V_i[iDim + idx.Velocity()] - GridVel_i[iDim];
          su2double Velocity_j = V_j[iDim + idx.Velocity()] - GridVel_j[iDim];
          q_ij += 0.5 * (Velocity_i + Velocity_j) * Normal[iDim];
        }
      } else {
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          q_ij += 0.5 * (V_i[iDim + idx.Velocity()] + V_j[iDim + idx.Velocity()]) * Normal[iDim];
        }
      }
      a0 = fmax(0.0, q_ij);
      a1 = fmin(0.0, q_ij);
    }

    FinishResidualCalc(config);

    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }
};
