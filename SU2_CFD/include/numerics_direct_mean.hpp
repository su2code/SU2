/*!
 * \file numerics_direct_mean.hpp
 * \brief
 * \author C. Pederson
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "numerics_structure.hpp"

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
  const unsigned short nPrimVar;  /*!< \brief The size of the primitive variable array used in the numerics class. */
  const bool correct_gradient; /*!< \brief Apply a correction to the gradient term */
  bool implicit;               /*!< \brief Implicit calculus. */
  su2double *heat_flux_vector, /*!< \brief Flux of total energy due to molecular and turbulent diffusion */
  *heat_flux_jac_i,            /*!< \brief Jacobian of the molecular + turbulent heat flux vector, projected onto the normal vector. */
  **tau_jacobian_i;            /*!< \brief Jacobian of the viscous + turbulent stress tensor, projected onto the normal vector. */
  su2double *Mean_PrimVar,     /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,      /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,          /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,                /*!< \brief Mean value of the turbulent kinetic energy. */
  Mean_TauWall,                /*!< \brief Mean wall shear stress (wall functions). */
  TauWall_i, TauWall_j,        /*!< \brief Wall shear stress at point i and j (wall functions). */
  dist_ij_2,                   /*!< \brief Length of the edge and face, squared */
  *Proj_Mean_GradPrimVar_Edge, /*!< \brief Inner product of the Mean gradient and the edge vector. */
  *Edge_Vector;                /*!< \brief Vector from point i to point j. */

  /*!
   * \brief Calculate the viscous + turbulent stress tensor
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  void GetStressTensor(const su2double *val_primvar,
                       su2double **val_gradprimvar,
                       const su2double val_turb_ke,
                       const su2double val_laminar_viscosity,
                       const su2double val_eddy_viscosity);

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
  void AddQCR(su2double **val_gradprimvar);

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
                  const su2double val_tau_wall);


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
  void GetTauJacobian(const su2double* val_Mean_PrimVar,
                      const su2double val_laminar_viscosity,
                      const su2double val_eddy_viscosity,
                      const su2double val_dist_ij,
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
  void GetIncTauJacobian(const su2double val_laminar_viscosity,
                         const su2double val_eddy_viscosity,
                         const su2double val_dist_ij,
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
                          const su2double val_dS,
                          const su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);

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
                bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Base();

  /*!
   * \brief Set the value of the wall shear stress at point i and j (wall functions).
   * \param[in] val_tauwall_i - Value of the wall shear stress at point i.
   * \param[in] val_tauwall_j - Value of the wall shear stress at point j.
   */
  void SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j);
};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 */
class CAvgGrad_Flow : public CAvgGrad_Base {
private:

  /*!
   * \brief Compute the heat flux due to molecular and turbulent diffusivity
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  void GetHeatFluxVector(su2double **val_gradprimvar,
                         const su2double val_laminar_viscosity,
                         const su2double val_eddy_viscosity);

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
  void GetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           const su2double val_laminar_viscosity,
                           const su2double val_eddy_viscosity,
                           const su2double val_dist_ij,
                           const su2double *val_normal);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CGeneralAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author M.Pini, S. Vitale
 */
class CGeneralAvgGrad_Flow : public CAvgGrad_Base {
private:
  su2double *Mean_SecVar,    /*!< \brief Mean secondary variables. */
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
  void GetHeatFluxVector(su2double **val_gradprimvar,
                          const su2double val_laminar_viscosity,
                          const su2double val_eddy_viscosity,
                          const su2double val_thermal_conductivity,
                          const su2double val_heat_capacity_cp);

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
  void GetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           const su2double *val_Mean_SecVar,
                           const su2double val_eddy_viscosity,
                           const su2double val_thermal_conductivity,
                           const su2double val_heat_capacity_cp,
                           const su2double val_dist_ij);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGrad_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradInc_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, F. Palacios, T. Economon
 */
class CAvgGradInc_Flow : public CAvgGrad_Base {
private:
  su2double Mean_Thermal_Conductivity, /*!< \brief Mean value of the effective thermal conductivity. */
  Mean_Mean_Cp,   /*!< \brief Mean value of the effective thermal conductivity and specific heat at constant pressure. */
  proj_vector_ij; /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  bool energy;    /*!< \brief computation with the energy equation. */

  /*
   * \brief Compute the projection of the viscous fluxes into a direction (artificial compresibility method).
   *
   * The viscous + turbulent stress tensor must be calculated before calling
   * this function.
   *
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   */
  void GetViscousIncProjFlux(su2double **val_gradprimvar,
                             su2double *val_normal,
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
                   bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradInc_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

// TODO: Move these to an inline file after refactoring at the end
inline void CAvgGrad_Base::SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) {
  TauWall_i = val_tauwall_i;
  TauWall_j = val_tauwall_j;
}
