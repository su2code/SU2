/*!
 * \file scalar_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in scalar problems.
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
 * \class CNoFlowIndices
 * \brief Dummy flow indices class to use CAvgGrad_Scalar when flow variables are not available.
 * For example, solid heat transfer problems.
 */
struct CNoFlowIndices {
  CNoFlowIndices(int, int) {}
  inline int Density() const { return 0; }
  inline int LaminarViscosity() const { return 0; }
  inline int EddyViscosity() const { return 0; }
};

/*!
 * \class CAvgGrad_Scalar
 * \brief Template class for computing viscous residual of scalar values
 * \details This class serves as a template for the scalar viscous residual
 *   classes.  The general structure of a viscous residual calculation is the
 *   same for many different  models, which leads to a lot of repeated code.
 *   By using the template design pattern, these sections of repeated code are
 *   moved to a shared base class, and the specifics of each model
 *   are implemented by derived classes.  In order to add a new residual
 *   calculation for a viscous residual, extend this class and implement
 *   the pure virtual functions with model-specific behavior.
 * \ingroup ViscDiscr
 * \author C. Pederson, A. Bueno, and F. Palacios
 */
template <class FlowIndices>
class CAvgGrad_Scalar : public CNumerics {
 protected:
  enum : unsigned short {MAXNVAR = 8};

  const FlowIndices idx;                      /*!< \brief Object to manage the access to the flow primitives. */
  su2double Proj_Mean_GradScalarVar[MAXNVAR]; /*!< \brief Mean_gradScalarVar DOT normal, corrected if required. */
  su2double proj_vector_ij = 0.0;             /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  su2double Flux[MAXNVAR];                    /*!< \brief Final result, diffusive flux/residual. */
  su2double* Jacobian_i[MAXNVAR];             /*!< \brief Flux Jacobian w.r.t. node i. */
  su2double* Jacobian_j[MAXNVAR];             /*!< \brief Flux Jacobian w.r.t. node j. */
  su2double JacobianBuffer[2*MAXNVAR*MAXNVAR];/*!< \brief Static storage for the two Jacobians. */

  const bool correct_gradient = false, incompressible = false;

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method, derived classes
   *        should compute the Flux and Jacobians (i/j) inside this method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(const CConfig* config) = 0;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_gradient - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Scalar(unsigned short val_nDim, unsigned short val_nVar, bool correct_grad,
                  const CConfig* config)
    : CNumerics(val_nDim, val_nVar, config),
      idx(val_nDim, config->GetnSpecies()),
      correct_gradient(correct_grad),
      incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) {
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
   * \brief Compute the viscous residual using an average of gradients without correction.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final {
    AD::StartPreacc();
    AD::SetPreaccIn(Coord_i, nDim);
    AD::SetPreaccIn(Coord_j, nDim);
    AD::SetPreaccIn(Normal, nDim);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(ScalarVar_Grad_j, nVar, nDim);
    if (correct_gradient) {
      AD::SetPreaccIn(ScalarVar_i, nVar);
      AD::SetPreaccIn(ScalarVar_j, nVar);
    }
    if (!std::is_same<FlowIndices, CNoFlowIndices>::value) {
      AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);
      AD::SetPreaccIn(V_j[idx.Density()], V_j[idx.LaminarViscosity()], V_j[idx.EddyViscosity()]);

      Density_i = V_i[idx.Density()];
      Density_j = V_j[idx.Density()];
      Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
      Laminar_Viscosity_j = V_j[idx.LaminarViscosity()];
      Eddy_Viscosity_i = V_i[idx.EddyViscosity()];
      Eddy_Viscosity_j = V_j[idx.EddyViscosity()];
    }

    ExtraADPreaccIn();

    su2double ProjGradScalarVarNoCorr[MAXNVAR];
    proj_vector_ij = ComputeProjectedGradient(nDim, nVar, Normal, Coord_i, Coord_j, ScalarVar_Grad_i, ScalarVar_Grad_j,
                                              correct_gradient, ScalarVar_i, ScalarVar_j, ProjGradScalarVarNoCorr,
                                              Proj_Mean_GradScalarVar);
    FinishResidualCalc(config);

    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }
};
