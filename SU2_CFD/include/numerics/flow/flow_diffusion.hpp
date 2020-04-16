/*!
 * \file flow_diffusion.hpp
 * \brief Delarations of numerics classes for viscous flux computation.
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
 * \class CAvgGrad_Base
 * \brief A base class for computing viscous terms using an average of gradients.
 * \details This is the base class for the numerics classes that compute the
 * viscous fluxes for the flow solvers (i.e. compressible or incompressible
 * Navier Stokes).  The actual numerics classes derive from this class.
 * This class is used to share functions and variables that are common to all
 * of the flow viscous numerics.  For example, the turbulent stress tensor
 * is computed identically for all three derived classes.
 * \ingroup ViscDiscr
 * \author C. Pederson, A. Bueno, F. Palacios, T. Economon
 */
class CAvgGrad_Base : public CNumerics {
protected:
  const unsigned short nPrimVar;          /*!< \brief The size of the primitive variable array used in the numerics class. */
  const bool correct_gradient;            /*!< \brief Apply a correction to the gradient term */
  bool implicit = false;                  /*!< \brief Implicit calculus. */
  su2double
  heat_flux_vector[MAXNDIM] = {0.0},      /*!< \brief Flux of total energy due to molecular and turbulent diffusion */
  *heat_flux_jac_i = nullptr,             /*!< \brief Jacobian of the molecular + turbulent heat flux vector, projected onto the normal vector. */
  **tau_jacobian_i = nullptr;             /*!< \brief Jacobian of the viscous + turbulent stress tensor, projected onto the normal vector. */
  su2double *Mean_PrimVar = nullptr;      /*!< \brief Mean primitive variables. */
  const su2double
  *PrimVar_i = nullptr,
  *PrimVar_j = nullptr;                   /*!< \brief Primitives variables at point i and j. */
  su2double **Mean_GradPrimVar = nullptr, /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                 /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                    /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,                           /*!< \brief Mean value of the turbulent kinetic energy. */
  Mean_TauWall,                           /*!< \brief Mean wall shear stress (wall functions). */
  TauWall_i, TauWall_j,                   /*!< \brief Wall shear stress at point i and j (wall functions). */
  dist_ij_2,                              /*!< \brief Length of the edge and face, squared */
  Edge_Vector[MAXNDIM] = {0.0},           /*!< \brief Vector from point i to point j. */
  *Proj_Mean_GradPrimVar_Edge = nullptr;  /*!< \brief Inner product of the Mean gradient and the edge vector. */
  bool TauWallFlag_i, TauWallFlag_j;      /*!< \brief Flag for Wall shear stress at point i and j (wall functions/models). */

  su2double** Jacobian_i = nullptr;       /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j = nullptr;       /*!< \brief The Jacobian w.r.t. point j after computation. */

  /*!
   * \brief Add a correction using a Quadratic Constitutive Relation
   *
   * This function requires that the stress tensor already be
   * computed using \ref GetStressTensor
   *
   * See: Spalart, P. R., "Strategies for Turbulence Modelling and
   * Simulation," International Journal of Heat and Fluid Flow, Vol. 21,
   * 2000, pp. 252-263
   *
   * \param[in] val_gradprimvar
   */
  void AddQCR(const su2double* const *val_gradprimvar);

  /*!
   * \brief Scale the stress tensor using a predefined wall stress.
   *
   * This function requires that the stress tensor already be
   * computed using \ref GetStressTensor
   *
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_tau_wall - The wall stress
   */
  void AddTauWall(const su2double *val_normal,
                  su2double val_tau_wall);

  /**
   * \brief Calculate the Jacobian of the viscous + turbulent stress tensor
   *
   * This function is intended only for the compressible flow solver.
   * This Jacobian is projected onto the normal vector, so it is of dimension
   * [nDim][nVar]
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void SetTauJacobian(const su2double* val_Mean_PrimVar,
                      su2double val_laminar_viscosity,
                      su2double val_eddy_viscosity,
                      su2double val_dist_ij,
                      const su2double *val_normal);


  /**
   * \brief Calculate the Jacobian of the viscous and turbulent stress tensor
   *
   * This function is intended only for the incompressible flow solver.
   * This Jacobian is projected onto the normal vector, so it is of dimension
   * [nDim][nVar]
   *
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void SetIncTauJacobian(su2double val_laminar_viscosity,
                         su2double val_eddy_viscosity,
                         su2double val_dist_ij,
                         const su2double *val_normal);

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   *
   * The heat flux vector and the stress tensor must be calculated before
   * calling this function.
   *
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void GetViscousProjFlux(const su2double *val_primvar,
                          const su2double *val_normal);

  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians.
   *
   * The Jacobians of the heat flux vector and the stress tensor must be
   * calculated before calling this function.
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                          su2double val_dS,
                          const su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);

  /*!
   * \brief Apply a correction to the gradient to reduce the truncation error
   *
   * \param[in] val_PrimVar_i - Primitive variables at point i
   * \param[in] val_PrimVar_j - Primitive variables at point j
   * \param[in] val_edge_vector - The vector between points i and j
   * \param[in] val_dist_ij_2 - The distance between points i and j, squared
   * \param[in] val_nPrimVar - The number of primitive variables
   */
  void CorrectGradient(su2double** GradPrimVar,
                       const su2double* val_PrimVar_i,
                       const su2double* val_PrimVar_j,
                       const su2double* val_edge_vector,
                       su2double val_dist_ij_2,
                       const unsigned short val_nPrimVar);

  /*!
   * \brief Initialize the Reynolds Stress Matrix
   * \param[in] turb_ke turbulent kinetic energy of node
   */
  void SetReynoldsStressMatrix(su2double turb_ke);

  /*!
   * \brief Perturb the Reynolds stress tensor based on parameters
   * \param[in] turb_ke: turbulent kinetic energy of the noce
   * \param[in] Eig_Val_Comp: Defines type of eigenspace perturbation
   * \param[in] beta_delta: Defines the amount of eigenvalue perturbation
   */
  void SetPerturbedRSM(su2double turb_ke, const CConfig* config);

  /*!
   * \brief Get the mean rate of strain matrix based on velocity gradients
   * \param[in] S_ij
   */
  void GetMeanRateOfStrainMatrix(su2double **S_ij) const;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables to use.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Base(unsigned short val_nDim, unsigned short val_nVar,
                unsigned short val_nPrimVar,
                bool val_correct_grad, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Base();

  /*!
   * \brief Set the value of the wall shear stress at point i and j (wall functions).
   * \param[in] val_tauwall_i - Value of the wall shear stress at point i.
   * \param[in] val_tauwall_j - Value of the wall shear stress at point j.
   */
  inline void SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) override {
    TauWall_i = val_tauwall_i;
    TauWall_j = val_tauwall_j;
  }


  /*!
   * \brief Set the value of the bollean flag to use (or not) the wall shear stress from the wall function.
   * \param[in] val_tauwallflag_i - Flag for Tauwall at point i
   * \param[in] val_tauwallflag_j - Flag for Tauwall at point j
   */
  inline void SetTauWall_Flag(bool val_tauwallflag_i, bool val_tauwallflag_j) override {
    TauWallFlag_i = val_tauwallflag_i;
    TauWallFlag_j = val_tauwallflag_j;
  }


  /*!
   * \brief Calculate the viscous + turbulent stress tensor
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  void SetStressTensor(const su2double *val_primvar,
                       const su2double* const *val_gradprimvar,
                       su2double val_turb_ke,
                       su2double val_laminar_viscosity,
                       su2double val_eddy_viscosity);

  /*!
   * \brief Get a component of the viscous stress tensor.
   *
   * \param[in] iDim - The first index
   * \param[in] jDim - The second index
   * \return The component of the viscous stress tensor at iDim, jDim
   */
  inline su2double GetStressTensor(unsigned short iDim, unsigned short jDim) const { return tau[iDim][jDim];}

  /*!
   * \brief Get a component of the heat flux vector.
   * \param[in] iDim - The index of the component
   * \return The component of the heat flux vector at iDim
   */
  inline su2double GetHeatFluxVector(unsigned short iDim) const { return heat_flux_vector[iDim]; }

};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 */
class CAvgGrad_Flow final : public CAvgGrad_Base {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar,
                bool val_correct_grad, const CConfig* config);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

  /*!
   * \brief Compute the heat flux due to molecular and turbulent diffusivity
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  void SetHeatFluxVector(const su2double* const *val_gradprimvar,
                         su2double val_laminar_viscosity,
                         su2double val_eddy_viscosity);

  /*!
   * \brief Compute the Jacobian of the heat flux vector
   *
   * This Jacobian is projected onto the normal vector, so it is of
   * dimension nVar.
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_gradprimvar - Mean value of the gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           su2double val_laminar_viscosity,
                           su2double val_eddy_viscosity,
                           su2double val_dist_ij,
                           const su2double *val_normal);
};

/*!
 * \class CAvgGradInc_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, F. Palacios, T. Economon
 */
class CAvgGradInc_Flow final : public CAvgGrad_Base {
private:
  su2double Mean_Thermal_Conductivity; /*!< \brief Mean value of the effective thermal conductivity. */
  bool energy;                         /*!< \brief computation with the energy equation. */

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction
   *
   * The viscous + turbulent stress tensor must be calculated before calling
   * this function.
   *
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   */
  void GetViscousIncProjFlux(const su2double* const *val_gradprimvar,
                             const su2double *val_normal,
                             su2double val_thermal_conductivity);

  /*!
   * \brief Compute the projection of the viscous Jacobian matrices.
   *
   * The Jacobian of the stress tensor must be calculated before calling
   * this function.
   *
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousIncProjJacs(su2double val_dS,
                             su2double **val_Proj_Jac_Tensor_i,
                             su2double **val_Proj_Jac_Tensor_j);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradInc_Flow(unsigned short val_nDim, unsigned short val_nVar,
                   bool val_correct_grad, const CConfig* config);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};


/*!
 * \class CGeneralAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author M.Pini, S. Vitale
 */
class CGeneralAvgGrad_Flow final : public CAvgGrad_Base {
private:
  su2double Mean_SecVar[2],  /*!< \brief Mean secondary variables. */
  Mean_Thermal_Conductivity, /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp;                   /*!< \brief Mean value of the Cp. */

  /*!
   * \brief Compute the heat flux due to molecular and turbulent diffusivity
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_heat_capacity_cp - Heat Capacity at constant pressure.
   */
  void SetHeatFluxVector(const su2double* const *val_gradprimvar,
                         su2double val_laminar_viscosity,
                         su2double val_eddy_viscosity,
                         su2double val_thermal_conductivity,
                         su2double val_heat_capacity_cp);

  /*!
   * \brief Compute the Jacobian of the heat flux vector
   *
   * This Jacobian is projected onto the normal vector, so it is of
   * dimension nVar.
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   */
  void SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           const su2double *val_Mean_SecVar,
                           su2double val_eddy_viscosity,
                           su2double val_thermal_conductivity,
                           su2double val_heat_capacity_cp,
                           su2double val_dist_ij);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, bool val_correct_grad, const CConfig* config);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};
