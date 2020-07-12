/*!
 * \file NEMO_diffusion.hpp
 * \brief Delarations of numerics classes for viscous flux computation. The implementation is in NEMO_diffusion.cpp.
 * \author F. Palacios, T. Economon
 * \version 7.0.5 "Blackbird"
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
#include "../../variables/CNEMOEulerVariable.hpp"

/*!
 * \class CAvgGrad_NEMO
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author S. Copeland, W. Maier
 * \version 6.2.0 "falcon"
 */
class CAvgGrad_NEMO : public CNumerics {
private:
  unsigned short iDim, iVar, nPrimVar, nPrimVarGrad;		/*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
  *Mean_U,
  **Mean_GU,
  *Mean_dTdU,
  *Mean_dTvedU,
  *Mean_dPdU,
  *Mean_Eve,
  *Mean_Cvve,
  *PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,						/*!< \brief Mean value of the gradient. */
  *Mean_Diffusion_Coeff, /*!< \brief Mean value of the species diffusion coefficient. */
  Mean_Laminar_Viscosity, /*!< \brief Mean value of the viscosity. */
  Mean_Thermal_Conductivity, /*!< \brief Mean value of the thermal conductivity. */
  Mean_Thermal_Conductivity_ve, /*!< \brief Mean value of the vib-el. thermal conductivity. */

  *ProjFlux,	/*!< \brief Projection of the viscous fluxes. */
  dist_ij;						/*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  CNEMOEulerVariable *variable;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of variables in the primitive variable gradient.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_NEMO(unsigned short val_nDim,
                unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_NEMO(void);

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of Primitive Variables.
   * \param[in] val_eve - Virbational-Electronical Energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_diffusioncoeff - Disffusion Coefficient.
   * \param[in] val_viscosity - Viscosity
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   * \param[in] val_thermal_conductivity_ve - Thermal conductivity of Vibe-Elec modes.
   * \param[in] config - Definition of the particular problem.
   */
  void GetViscousProjFlux(su2double *val_primvar,
                          su2double **val_gradprimvar,
                          su2double *val_eve,
                          const su2double *val_normal,
                          su2double *val_diffusioncoeff,
                          su2double val_viscosity,
                          su2double val_therm_conductivity,
                          su2double val_therm_conductivity_ve,
                          CConfig *config);
  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians for arbitrary equations of state.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_gradprimvar - Mean value of the gradient of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(su2double *val_Mean_PrimVar,
                          su2double **val_Mean_GradPrimVar,
                          su2double *val_Mean_Eve,
                          su2double *val_Mean_Cvve,
                          su2double *val_diffusion_coeff,
                          su2double val_laminar_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_thermal_conductivity_ve,
                          su2double val_dist_ij, su2double *val_normal,
                          su2double val_dS, su2double *val_Fv,
                          su2double **val_Jac_i, su2double **val_Jac_j,
                          CConfig *config);
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i,
                       su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CAvgGradCorrected_NEMO
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author S. Copeland, W. Maier
 * \version 6.1.0 "Falcon"
 */
class CAvgGradCorrected_NEMO : public CNumerics {
private:
  unsigned short iDim, iVar, nPrimVar, nPrimVarGrad;		/*!< \brief Iterators in dimension an variable. */
  su2double
  *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,						/*!< \brief Mean value of the gradient. */
  *Mean_Eve,
  *Mean_Cvve,
  *Edge_Vector,
  *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  *Mean_Diffusion_Coeff, /*!< \brief Mean value of the species diffusion coefficient. */
  Mean_Laminar_Viscosity, /*!< \brief Mean value of the viscosity. */
  Mean_Thermal_Conductivity, /*!< \brief Mean value of the thermal conductivity. */
  Mean_Thermal_Conductivity_ve, /*!< \brief Mean value of the vib-el. thermal conductivity. */

  *ProjFlux,	/*!< \brief Projection of the viscous fluxes. */
  dist_ij;						/*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */

  CNEMOEulerVariable *variable;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of variables in the primitive variable gradient.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_NEMO(unsigned short val_nDim,
                unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_NEMO(void);

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of Primitive Variables.
   * \param[in] val_eve - Virbational-Electronical Energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_diffusioncoeff - Disffusion Coefficient.
   * \param[in] val_viscosity - Viscosity
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   * \param[in] val_thermal_conductivity_ve - Thermal conductivity of Vibe-Elec modes.
   * \param[in] config - Definition of the particular problem.
   */
  void GetViscousProjFlux(su2double *val_primvar,
                          su2double **val_gradprimvar,
                          su2double *val_eve,
                          const su2double *val_normal,
                          su2double *val_diffusioncoeff,
                          su2double val_viscosity,
                          su2double val_therm_conductivity,
                          su2double val_therm_conductivity_ve,
                          CConfig *config);

  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians for arbitrary equations of state.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_gradprimvar - Mean value of the gradient of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(su2double *val_Mean_PrimVar,
                          su2double **val_Mean_GradPrimVar,
                          su2double *val_Mean_Eve,
                          su2double *val_Mean_Cvve,
                          su2double *val_diffusion_coeff,
                          su2double val_laminar_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_thermal_conductivity_ve,
                          su2double val_dist_ij, su2double *val_normal,
                          su2double val_dS, su2double *val_Fv,
                          su2double **val_Jac_i, su2double **val_Jac_j,
                          CConfig *config);
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i,
                       su2double **val_Jacobian_j,
                       CConfig *config);
};
