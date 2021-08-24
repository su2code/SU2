/*!
 * \file CNEMONumerics.hpp
 * \brief Base class template NEMO numerics.
 * \author C. Garbacz, W. Maier, S. R. Copeland
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CNEMONumerics
 * \brief Base class template NEMO numerics.
 * \author C. Garbacz.
 */
class CNEMONumerics : public CNumerics {
public:
  bool implicit, ionization;
  su2double *rhos_i, *u_i;
  su2double *rhos_j, *u_j;
  su2double a_i, P_i, h_i;
  su2double a_j, P_j, h_j;
  unsigned short nPrimVar, nPrimVarGrad;

  su2double* Flux = nullptr;            /*!< \brief The flux / residual across the edge. */

  unsigned short nSpecies, nHeavy, nEl; /*!< \brief Number of species present in plasma */

  su2double *dPdU_i, *dPdU_j;
  su2double *dTdU_i, *dTdU_j;
  su2double *dTvedU_i, *dTvedU_j;
  su2double Gamma_i, Gamma_j;

  vector<su2double> hs;
  su2double *eve_i, *eve_j, *Cvve_i, *Cvve_j;

  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX,
  LAM_VISC_INDEX, EDDY_VISC_INDEX;

  CNEMOGas *fluidmodel;

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of primitive grad. variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMONumerics(unsigned short val_nDim, unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                const CConfig* config);

  /*!
   * \brief Constructor of the NEMO class, for turbulence applications.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMONumerics(unsigned short val_nDim, unsigned short val_nVar,
                unsigned short val_nVar_NEMO, unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CNEMONumerics(void);

  /*!
   * \Overload
   * \brief Compute the projected inviscid flux vector.
   * \param[in] val_U - Pointer to the conserved variables.
   * \param[in] val_V - Pointer to the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidProjFlux(const su2double *val_U, const su2double *val_V,
                           const su2double *val_normal, su2double *val_Proj_Flux);

  /*!
   * \overload
   * \brief Compute the projection of the inviscid Jacobian matrices for the two-temperature model.
   * \param[in] val_U - Vector conserved variables.
   * \param[in] val_V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of partial derivatives of pressure w.r.t. conserved vars.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(const su2double *val_U, const su2double *val_V, const su2double *val_dPdU,
                          const su2double *val_normal, const su2double val_scale,
                          su2double **val_Proj_Jac_Tensor);

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of Primitive Variables.
   * \param[in] val_eve - Virbational-Electronical Energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_diffusioncoeff - Disffusion Coefficient.
   * \param[in] val_lam_viscosity - Laminar Viscosity
   * \param[in] val_eddy_viscosity - Eddy Viscosity
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   * \param[in] val_thermal_conductivity_ve - Thermal conductivity of Vibe-Elec modes.
   * \param[in] config - Definition of the particular problem.
   */
  void GetViscousProjFlux(su2double *val_primvar,
                          su2double **val_gradprimvar,
                          su2double *val_eve,
                          const su2double *val_normal,
                          su2double *val_diffusioncoeff,
                          su2double val_lam_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_therm_conductivity,
                          su2double val_therm_conductivity_ve,
                          const CConfig *config);
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
                          su2double val_eddy_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_thermal_conductivity_ve,
                          su2double val_dist_ij, su2double *val_normal,
                          su2double val_dS, su2double *val_Fv,
                          su2double **val_Jac_i, su2double **val_Jac_j,
                          const CConfig *config);

  /*!
   * \overload
   * \brief Computation of the matrix P, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables (really only need rhoEve)
   * \param[in] V - Vector of primitive variables
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved vars.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix(const su2double *U, const su2double *V, const su2double *val_dPdU,
                  const su2double *val_normal, const su2double *l, const su2double *m,
                  su2double **val_p_tensor) const;

  /*!
   * \overload
   * \brief Computation of the matrix P^{-1}, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables.
   * \param[in] V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved variables
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(const su2double *U, const su2double *V, const su2double *val_dPdU,
                      const su2double *val_normal, const su2double *l, const su2double *m,
                      su2double **val_invp_tensor) const;


  /*!
   * \brief Set the pressure derivatives.
   * \param[in] val_dPdU_i - pressure derivatives at i.
   * \param[in] val_dPdU_j - pressure derivatives at j.
   */
  inline void SetdPdU(su2double *val_dPdU_i, su2double *val_dPdU_j)       final { dPdU_i = val_dPdU_i; dPdU_j = val_dPdU_j; }

  /*!
   * \brief Set the temperature derivatives.
   * \param[in] val_dTdU_i - temperature derivatives at i.
   * \param[in] val_dTdU_j - temperature derivatives at j.
   */
  inline void SetdTdU(su2double *val_dTdU_i, su2double *val_dTdU_j)       final { dTdU_i = val_dTdU_i; dTdU_j = val_dTdU_j; }

  /*!
   * \brief Set the vib-el temperature derivatives.
   * \param[in] val_dTvedU_i - t_ve derivatives at i.
   * \param[in] val_dTvedU_j - t_ve derivatives at j.
   */
  inline void SetdTvedU(su2double *val_dTvedU_i, su2double *val_dTvedU_j) final { dTvedU_i = val_dTvedU_i; dTvedU_j = val_dTvedU_j; }

  /*!
   * \brief Set the vib-el energy.
   * \param[in] val_Eve_i - vib-el energy at i.
   * \param[in] val_Eve_j - vib-el energy at j.
   */
  inline void SetEve(su2double *val_Eve_i, su2double *val_Eve_j)          final {eve_i = val_Eve_i; eve_j = val_Eve_j; }

  /*!
   * \brief Set the Cvve.
   * \param[in] val_Cvve_i - cvve at i.
   * \param[in] val_Cvve_j - cvve at j.
   */
  inline void SetCvve(su2double *val_Cvve_i, su2double *val_Cvve_j)       final {Cvve_i = val_Cvve_i; Cvve_j = val_Cvve_j; }

  /*!
   * \brief Set the ratio of specific heats.
   * \param[in] val_Gamma_i - Gamma at i.
   * \param[in] val_Gamma_j - Gamma at j.
   */
  inline void SetGamma(su2double val_Gamma_i, su2double val_Gamma_j)      final {Gamma_i = val_Gamma_i; Gamma_j = val_Gamma_j; }

};
