/*!
 * \file turb_convection.hpp
 * \brief Delarations of numerics classes for discretization of
 *        convective fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
class CUpwScalar : public CNumerics {
protected:
  su2double
  Lambda[3] = {0.0},           /*!< \brief Projected velocity at the face. */
  Epsilon[3] = {0.0},          /*!< \brief Band for entropy correction. */
  ProjVel_i = 0.0,             /*!< \brief The maximum of the face-normal velocity and 0 */
  ProjVel_j = 0.0,             /*!< \brief The minimum of the face-normal velocity and 0 */
  SoundSpeed_i = 0.0,
  SoundSpeed_j = 0.0,
  *Flux = nullptr,             /*!< \brief Final result, diffusive flux/residual. */
  **Jacobian_i = nullptr,      /*!< \brief Flux Jacobian w.r.t. node i. */
  **Jacobian_j = nullptr;      /*!< \brief Flux Jacobian w.r.t. node j. */
  su2double muscl_kappa;
  bool muscl;

  su2double Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
            Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j;
  su2double Velocity_i[MAXNDIM] = {0.0}, Velocity_j[MAXNDIM] = {0.0};
  su2double R, R_Plus_One, RoeSqVel, RoeEnthalpy, RoeTke, RoeOmega, RoeSoundSpeed2, RoeSoundSpeed;

  const bool implicit = false, incompressible = false, dynamic_grid = false;

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Compute the contribution of central differencing to the flux Jacobian.
   * \param[in] val_kappa - MUSCL kappa blending parameter.
   * \param[in/out] jac_{i/j} - Flux Jacobian wrt conservative at [face/node].
   * \param[in] lim_{i/j} - Slope limiter at node.
   * \param[in] r_{i/j} - Density at face.
   * \param[in] r_n_{i/j} - Density at node.
   */
  void GetMUSCLJac(su2double **jac_i, su2double **jac_j,
                   const su2double *lim_i, const su2double *lim_j,
                   const su2double *r_i, const su2double *r_j,
                   const su2double *r_n_i, const su2double *r_n_j);

  /*!
   * \brief Model-specific steps in the ComputeResidual method, derived classes
   *        compute the Flux and its Jacobians via this method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(const CConfig* config) = 0;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwScalar(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_muscl);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwScalar(void) override;

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaras turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 */
class CUpwSca_TurbSA final : public CUpwScalar {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_muscl);

};

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 */
class CUpwSca_TurbSST final : public CUpwScalar {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_muscl);

};
