/*!
 * \file CNumerics.hpp
 * \brief Delaration of the base numerics class, the
 *        implementation is in the CNumerics.cpp file.
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

#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../../Common/include/CConfig.hpp"

using namespace std;

class CElement;
class CFluidModel;

/*!
 * \class CNumerics
 * \brief Class for defining the numerical methods.
 * \author F. Palacios, T. Economon
 */
class CNumerics {
protected:
  enum : size_t {MAXNDIM = 3}; /*!< \brief Max number of space dimensions, used in some static arrays. */

  unsigned short nDim, nVar;  /*!< \brief Number of dimensions and variables. */
  su2double Gamma;            /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;  /*!< \brief Fluids's Gamma - 1.0  . */
  su2double Minf;             /*!< \brief Free stream Mach number . */
  su2double Gas_Constant;     /*!< \brief Gas constant. */
  su2double *Vector;          /*!< \brief Auxiliary vector. */
  su2double Prandtl_Lam;      /*!< \brief Laminar Prandtl's number. */
  su2double Prandtl_Turb;     /*!< \brief Turbulent Prandtl's number. */
  su2double
  **Flux_Tensor,      /*!< \brief Flux tensor (used for viscous and inviscid purposes. */
  *Proj_Flux_Tensor;  /*!< \brief Flux tensor projected in a direction. */
  su2double
  **tau,      /*!< \brief Viscous stress tensor. */
  **delta,    /*!< \brief Identity matrix. */
  **delta3;   /*!< \brief 3 row Identity matrix. */
  su2double
  *Diffusion_Coeff_i, /*!< \brief Species diffusion coefficients at point i. */
  *Diffusion_Coeff_j; /*!< \brief Species diffusion coefficients at point j. */
  su2double
  Laminar_Viscosity_i,   /*!< \brief Laminar viscosity at point i. */
  Laminar_Viscosity_j,   /*!< \brief Laminar viscosity at point j. */
  Laminar_Viscosity_id,  /*!< \brief Variation of laminar viscosity at point i. */
  Laminar_Viscosity_jd;  /*!< \brief Variation of laminar viscosity at point j. */
  su2double
  Thermal_Conductivity_i,    /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_j,    /*!< \brief Thermal conductivity at point j. */
  Thermal_Diffusivity_i,     /*!< \brief Thermal diffusivity at point i. */
  Thermal_Diffusivity_j;     /*!< \brief Thermal diffusivity at point j. */
  su2double
  Cp_i,               /*!< \brief Cp at point i. */
  Cp_j;               /*!< \brief Cp at point j. */
  su2double
  Eddy_Viscosity_i,  /*!< \brief Eddy viscosity at point i. */
  Eddy_Viscosity_j;  /*!< \brief Eddy viscosity at point j. */
  su2double
  turb_ke_i,  /*!< \brief Turbulent kinetic energy at point i. */
  turb_ke_j;  /*!< \brief Turbulent kinetic energy at point j. */
  su2double
  Pressure_i,  /*!< \brief Pressure at point i. */
  Pressure_j;  /*!< \brief Pressure at point j. */
  su2double
  GravityForce_i,  /*!< \brief Gravity force at point i. */
  GravityForce_j;  /*!< \brief Gravity force at point j. */
  su2double
  Density_i,  /*!< \brief Density at point i. */
  Density_j;  /*!< \brief Density at point j. */
  su2double
  DensityInc_i,  /*!< \brief Incompressible density at point i. */
  DensityInc_j;  /*!< \brief Incompressible density at point j. */
  su2double
  BetaInc2_i,  /*!< \brief Beta incompressible at point i. */
  BetaInc2_j;  /*!< \brief Beta incompressible at point j. */
  su2double
  Lambda_i,  /*!< \brief Spectral radius at point i. */
  Lambda_j;  /*!< \brief Spectral radius at point j. */
  su2double
  LambdaComb_i,  /*!< \brief Spectral radius at point i. */
  LambdaComb_j;  /*!< \brief Spectral radius at point j. */
  su2double
  SoundSpeed_i,  /*!< \brief Sound speed at point i. */
  SoundSpeed_j;  /*!< \brief Sound speed at point j. */
  su2double
  Enthalpy_i,  /*!< \brief Enthalpy at point i. */
  Enthalpy_j;  /*!< \brief Enthalpy at point j. */
  su2double
  dist_i,  /*!< \brief Distance of point i to the nearest wall. */
  dist_j;  /*!< \brief Distance of point j to the nearest wall. */
  su2double
  Temp_i,  /*!< \brief Temperature at point i. */
  Temp_j;  /*!< \brief Temperature at point j. */
  su2double
  *Und_Lapl_i,  /*!< \brief Undivided laplacians at point i. */
  *Und_Lapl_j;  /*!< \brief Undivided laplacians at point j. */
  su2double
  Sensor_i,  /*!< \brief Pressure sensor at point i. */
  Sensor_j;  /*!< \brief Pressure sensor at point j. */
  su2double
  *GridVel_i,  /*!< \brief Grid velocity at point i. */
  *GridVel_j;  /*!< \brief Grid velocity at point j. */
  su2double
  *U_i,           /*!< \brief Vector of conservative variables at point i. */
  *U_id,          /*!< \brief Vector of derivative of conservative variables at point i. */
  *U_j,           /*!< \brief Vector of conservative variables at point j. */
  *U_jd;          /*!< \brief Vector of derivative of conservative variables at point j. */
  su2double
  *V_i,     /*!< \brief Vector of primitive variables at point i. */
  *V_j;     /*!< \brief Vector of primitive variables at point j. */
  su2double
  *S_i,     /*!< \brief Vector of secondary variables at point i. */
  *S_j;     /*!< \brief Vector of secondary variables at point j. */
  su2double
  *Psi_i,    /*!< \brief Vector of adjoint variables at point i. */
  *Psi_j;    /*!< \brief Vector of adjoint variables at point j. */
  su2double
  *DeltaU_i,  /*!< \brief Vector of linearized variables at point i. */
  *DeltaU_j;  /*!< \brief Vector of linearized variables at point j. */
  su2double
  *TurbVar_i,   /*!< \brief Vector of turbulent variables at point i. */
  *TurbVar_id,  /*!< \brief Vector of derivative of turbulent variables at point i. */
  *TurbVar_j,   /*!< \brief Vector of turbulent variables at point j. */
  *TurbVar_jd;  /*!< \brief Vector of derivative of turbulent variables at point j. */
  su2double
  *TransVar_i,  /*!< \brief Vector of turbulent variables at point i. */
  *TransVar_j;  /*!< \brief Vector of turbulent variables at point j. */
  su2double
  *TurbPsi_i,  /*!< \brief Vector of adjoint turbulent variables at point i. */
  *TurbPsi_j;  /*!< \brief Vector of adjoint turbulent variables at point j. */
  su2double
  **ConsVar_Grad_i,  /*!< \brief Gradient of conservative variables at point i. */
  **ConsVar_Grad_j,  /*!< \brief Gradient of conservative variables at point j. */
  **ConsVar_Grad;    /*!< \brief Gradient of conservative variables which is a scalar. */
  su2double
  **PrimVar_Grad_i,  /*!< \brief Gradient of primitive variables at point i. */
  **PrimVar_Grad_j;  /*!< \brief Gradient of primitive variables at point j. */
  su2double
  **PsiVar_Grad_i,  /*!< \brief Gradient of adjoint variables at point i. */
  **PsiVar_Grad_j;  /*!< \brief Gradient of adjoint variables at point j. */
  su2double
  **TurbVar_Grad_i,  /*!< \brief Gradient of turbulent variables at point i. */
  **TurbVar_Grad_j;  /*!< \brief Gradient of turbulent variables at point j. */
  su2double
  **TransVar_Grad_i,  /*!< \brief Gradient of turbulent variables at point i. */
  **TransVar_Grad_j;  /*!< \brief Gradient of turbulent variables at point j. */
  su2double
  **TurbPsi_Grad_i,  /*!< \brief Gradient of adjoint turbulent variables at point i. */
  **TurbPsi_Grad_j;  /*!< \brief Gradient of adjoint turbulent variables at point j. */
  su2double
  *AuxVar_Grad_i,    /*!< \brief Gradient of an auxiliary variable at point i. */
  *AuxVar_Grad_j;    /*!< \brief Gradient of an auxiliary variable at point i. */
  const su2double *RadVar_Source;  /*!< \brief Source term from the radiative heat transfer equation. */
  su2double
  *Coord_i,      /*!< \brief Cartesians coordinates of point i. */
  *Coord_j;      /*!< \brief Cartesians coordinates of point j. */
  unsigned short
  Neighbor_i,  /*!< \brief Number of neighbors of the point i. */
  Neighbor_j;  /*!< \brief Number of neighbors of the point j. */
  su2double
  *Normal,       /*!< \brief Normal vector, its norm is the area of the face. */
  *UnitNormal,   /*!< \brief Unitary normal vector. */
  *UnitNormald;  /*!< \brief Derivative of unitary normal vector. */
  su2double
  TimeStep,    /*!< \brief Time step useful in dual time method. */
  Area,        /*!< \brief Area of the face i-j. */
  Volume;      /*!< \brief Volume of the control volume around point i. */
  su2double vel2_inf;     /*!< \brief value of the square of freestream speed. */
  su2double
  *WindGust_i,  /*!< \brief Wind gust at point i. */
  *WindGust_j;  /*!< \brief Wind gust at point j. */
  su2double
  *WindGustDer_i,  /*!< \brief Wind gust derivatives at point i. */
  *WindGustDer_j;  /*!< \brief Wind gust derivatives at point j. */
  su2double *Vorticity_i, *Vorticity_j;    /*!< \brief Vorticity. */
  su2double StrainMag_i, StrainMag_j;      /*!< \brief Strain rate magnitude. */
  su2double Dissipation_i, Dissipation_j;  /*!< \brief Dissipation. */
  su2double Dissipation_ij;

  su2double *l, *m;

  su2double **MeanReynoldsStress; /*!< \brief Mean Reynolds stress tensor  */
  su2double **MeanPerturbedRSM;   /*!< \brief Perturbed Reynolds stress tensor  */
  bool using_uq;                  /*!< \brief Flag for UQ methodology  */
  su2double PerturbedStrainMag;   /*!< \brief Strain magnitude calculated using perturbed stress tensor  */
  unsigned short Eig_Val_Comp;    /*!< \brief Component towards which perturbation is perfromed */
  su2double uq_delta_b;           /*!< \brief Magnitude of perturbation */
  su2double uq_urlx;              /*!< \brief Under-relaxation factor for numerical stability */
  bool uq_permute;                /*!< \brief Flag for eigenvector permutation */

  /* Supporting data structures for the eigenspace perturbation for UQ methodology */
  su2double **A_ij, **newA_ij, **Eig_Vec, **New_Eig_Vec, **Corners;
  su2double *Eig_Val, *Barycentric_Coord, *New_Coord;

public:
  /*!
   * \brief Return type used in some "ComputeResidual" overloads to give a
   * const-view of the internally stored flux vector and Jacobians to the outside.
   * \note The default template params make this an "all const" object, it cannot
   * change after instantiation, nor can it be used to change the data.
   */
  template<class Vector_t = const su2double*,
           class Matrix_t = const Vector_t*>
  struct ResidualType {
    const Vector_t residual;
    const Matrix_t jacobian_i;
    const Matrix_t jacobian_j;

    ResidualType() = delete;

    ResidualType(const Vector_t& res, const Matrix_t& jac_i, const Matrix_t& jac_j) :
      residual(res), jacobian_i(jac_i), jacobian_j(jac_j) { }

    /*!
     * \brief The object can be directly cast to the vector type, this
     * allows discarding the Jacobians when they are not needed.
     */
    operator Vector_t() { return residual; }
  };

  /*!
   * \brief Constructor of the class.
   */
  CNumerics(void);

  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNumerics(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CNumerics(void);

  /*!
   * \brief Compute the determinant of a 3 by 3 matrix.
   * \param[in] val_matrix 3 by 3 matrix.
   * \return Determinant of the matrix
   */
  inline static su2double Determinant_3x3(su2double A00, su2double A01, su2double A02,
                                          su2double A10, su2double A11, su2double A12,
                                          su2double A20, su2double A21, su2double A22) {
    return A00*(A11*A22-A12*A21) - A01*(A10*A22-A12*A20) + A02*(A10*A21-A11*A20);
  }

  /*!
   * \brief Set the time step.
   * \param[in] val_timestep - Value of the time step.
   */
  inline void SetTimeStep(su2double val_timestep) { TimeStep = val_timestep;}

  /*!
   * \brief Set the freestream velocity square.
   * \param[in] SetVelocity2_Inf - Value of the square of the freestream velocity.
   */
  inline void SetVelocity2_Inf(su2double val_velocity2) { vel2_inf = val_velocity2; }

  /*!
   * \brief Set the value of the vorticity
   * \param[in] val_vorticity - Value of the vorticity.
   */
  void SetVorticity(su2double *val_vorticity_i, su2double *val_vorticity_j) {
    Vorticity_i = val_vorticity_i;
    Vorticity_j = val_vorticity_j;
  }

  /*!
   * \brief Set the value of the rate of strain magnitude.
   * \param[in] val_StrainMag_i - Value of the magnitude of rate of strain at point i.
   * \param[in] val_StrainMag_j - Value of the magnitude of rate of strain at point j.
   */
  void SetStrainMag(su2double val_strainmag_i, su2double val_strainmag_j) {
    StrainMag_i = val_strainmag_i;
    StrainMag_j = val_strainmag_j;
  }

  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_i - Value of the conservative variable at point i.
   * \param[in] val_u_j - Value of the conservative variable at point j.
   */
  inline void SetConservative(su2double *val_u_i, su2double *val_u_j) {
    U_i = val_u_i;
    U_j = val_u_j;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  inline void SetPrimitive(su2double *val_v_i, su2double *val_v_j) {
    V_i = val_v_i;
    V_j = val_v_j;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  inline void SetSecondary(su2double *val_s_i, su2double *val_s_j) {
    S_i = val_s_i;
    S_j = val_s_j;
  }

  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_i - Gradient of the conservative variable at point i.
   * \param[in] val_consvar_grad_j - Gradient of the conservative variable at point j.
   */
  inline void SetConsVarGradient(su2double **val_consvar_grad_i,
                                 su2double **val_consvar_grad_j) {
    ConsVar_Grad_i = val_consvar_grad_i;
    ConsVar_Grad_j = val_consvar_grad_j;
  }

  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad - Gradient of the conservative variable which is a scalar.
   */
  inline void SetConsVarGradient(su2double **val_consvar_grad) { ConsVar_Grad = val_consvar_grad; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
   * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
   */
  void SetPrimVarGradient(su2double **val_primvar_grad_i,
                          su2double **val_primvar_grad_j) {
    PrimVar_Grad_i = val_primvar_grad_i;
    PrimVar_Grad_j = val_primvar_grad_j;
  }

  /*!
   * \brief Set the value of the adjoint variable.
   * \param[in] val_psi_i - Value of the adjoint variable at point i.
   * \param[in] val_psi_j - Value of the adjoint variable at point j.
   */
  inline void SetAdjointVar(su2double *val_psi_i, su2double *val_psi_j) {
    Psi_i = val_psi_i;
    Psi_j = val_psi_j;
  }

  /*!
   * \brief Set the gradient of the adjoint variables.
   * \param[in] val_psivar_grad_i - Gradient of the adjoint variable at point i.
   * \param[in] val_psivar_grad_j - Gradient of the adjoint variable at point j.
   */
  inline void SetAdjointVarGradient(su2double **val_psivar_grad_i,
                                    su2double **val_psivar_grad_j) {
    PsiVar_Grad_i = val_psivar_grad_i;
    PsiVar_Grad_j = val_psivar_grad_j;
  }

  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
   * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
   */
  inline void SetTurbVar(su2double *val_turbvar_i, su2double *val_turbvar_j) {
    TurbVar_i = val_turbvar_i;
    TurbVar_j = val_turbvar_j;
  }

  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_transvar_i - Value of the turbulent variable at point i.
   * \param[in] val_transvar_j - Value of the turbulent variable at point j.
   */
  inline void SetTransVar(su2double *val_transvar_i, su2double *val_transvar_j) {
    TransVar_i = val_transvar_i;
    TransVar_j = val_transvar_j;
  }

  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  inline void SetTurbVarGradient(su2double **val_turbvar_grad_i,
                                 su2double **val_turbvar_grad_j) {
    TurbVar_Grad_i = val_turbvar_grad_i;
    TurbVar_Grad_j = val_turbvar_grad_j;
  }

  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  inline void SetTransVarGradient(su2double **val_transvar_grad_i,
                                  su2double **val_transvar_grad_j) {
    TransVar_Grad_i = val_transvar_grad_i;
    TransVar_Grad_j = val_transvar_grad_j;
  }

  /*!
   * \brief Set the value of the adjoint turbulent variable.
   * \param[in] val_turbpsivar_i - Value of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_j - Value of the adjoint turbulent variable at point j.
   */
  inline void SetTurbAdjointVar(su2double *val_turbpsivar_i, su2double *val_turbpsivar_j) {
    TurbPsi_i = val_turbpsivar_i;
    TurbPsi_j = val_turbpsivar_j;
  }

  /*!
   * \brief Set the gradient of the adjoint turbulent variables.
   * \param[in] val_turbpsivar_grad_i - Gradient of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_grad_j - Gradient of the adjoint turbulent variable at point j.
   */
  inline void SetTurbAdjointGradient(su2double **val_turbpsivar_grad_i,
                                     su2double **val_turbpsivar_grad_j) {
    TurbPsi_Grad_i = val_turbpsivar_grad_i;
    TurbPsi_Grad_j = val_turbpsivar_grad_j;
  }

  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first Menter blending function at point i.
   * \param[in] val_F1_j - Value of the first Menter blending function at point j.
   */
  virtual void SetF1blending(su2double val_F1_i, su2double val_F1_j) {/* empty */};

  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F1_i - Value of the second Menter blending function at point i.
   * \param[in] val_F1_j - Value of the second Menter blending function at point j.
   */
  virtual void SetF2blending(su2double val_F1_i, su2double val_F1_j) {/* empty */};

  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
   */
  virtual void SetCrossDiff(su2double val_CDkw_i, su2double val_CDkw_j) {/* empty */};

  /*!
   * \brief Set the gradient of the auxiliary variables.
   * \param[in] val_auxvargrad_i - Gradient of the auxiliary variable at point i.
   * \param[in] val_auxvargrad_j - Gradient of the auxiliary variable at point j.
   */
  inline void SetAuxVarGrad(su2double *val_auxvargrad_i, su2double *val_auxvargrad_j) {
    AuxVar_Grad_i = val_auxvargrad_i;
    AuxVar_Grad_j = val_auxvargrad_j;
  }

  /*!
   * \brief Set the diffusion coefficient
   * \param[in] val_diffusioncoeff_i - Value of the diffusion coefficients at i.
   * \param[in] val_diffusioncoeff_j - Value of the diffusion coefficients at j
   */
  inline void SetDiffusionCoeff(su2double* val_diffusioncoeff_i,
                                su2double* val_diffusioncoeff_j) {
    Diffusion_Coeff_i = val_diffusioncoeff_i;
    Diffusion_Coeff_j = val_diffusioncoeff_j;
  }

  /*!
   * \brief Set the laminar viscosity.
   * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
   * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
   */
  inline void SetLaminarViscosity(su2double val_laminar_viscosity_i,
                                  su2double val_laminar_viscosity_j) {
    Laminar_Viscosity_i = val_laminar_viscosity_i;
    Laminar_Viscosity_j = val_laminar_viscosity_j;
  }

  /*!
   * \brief Set the thermal conductivity (translational/rotational)
   * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
   * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  inline void SetThermalConductivity(su2double val_thermal_conductivity_i,
                                     su2double val_thermal_conductivity_j) {
    Thermal_Conductivity_i = val_thermal_conductivity_i;
    Thermal_Conductivity_j = val_thermal_conductivity_j;
  }

  /*!
   * \brief Set the thermal diffusivity (translational/rotational)
   * \param[in] val_thermal_diffusivity_i - Value of the thermal diffusivity at point i.
   * \param[in] val_thermal_diffusivity_j - Value of the thermal diffusivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  inline void SetThermalDiffusivity(su2double val_thermal_diffusivity_i,
                                    su2double val_thermal_diffusivity_j) {
    Thermal_Diffusivity_i = val_thermal_diffusivity_i;
    Thermal_Diffusivity_j = val_thermal_diffusivity_j;
  }

  /*!
   * \brief Set the eddy viscosity.
   * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
   * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
   */
  inline void SetEddyViscosity(su2double val_eddy_viscosity_i,
                               su2double val_eddy_viscosity_j) {
    Eddy_Viscosity_i = val_eddy_viscosity_i;
    Eddy_Viscosity_j = val_eddy_viscosity_j;
  }

  /*!
   * \brief Set the turbulent kinetic energy.
   * \param[in] val_turb_ke_i - Value of the turbulent kinetic energy at point i.
   * \param[in] val_turb_ke_j - Value of the turbulent kinetic energy at point j.
   */
  inline void SetTurbKineticEnergy(su2double val_turb_ke_i, su2double val_turb_ke_j) {
    turb_ke_i = val_turb_ke_i;
    turb_ke_j = val_turb_ke_j;
  }

  /*!
   * \brief Set the value of the distance from the nearest wall.
   * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
   * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
   */
  void SetDistance(su2double val_dist_i, su2double val_dist_j) {
    dist_i = val_dist_i;
    dist_j = val_dist_j;
  }

  /*!
   * \brief Set coordinates of the points.
   * \param[in] val_coord_i - Coordinates of the point i.
   * \param[in] val_coord_j - Coordinates of the point j.
   */
  inline void SetCoord(su2double *val_coord_i, su2double *val_coord_j) {
    Coord_i = val_coord_i;
    Coord_j = val_coord_j;
  }

  /*!
   * \brief Set the velocity of the computational grid.
   * \param[in] val_gridvel_i - Grid velocity of the point i.
   * \param[in] val_gridvel_j - Grid velocity of the point j.
   */
  inline void SetGridVel(su2double *val_gridvel_i, su2double *val_gridvel_j) {
    GridVel_i = val_gridvel_i;
    GridVel_j = val_gridvel_j;
  }

  /*!
   * \brief Set the wind gust value.
   * \param[in] val_windgust_i - Wind gust of the point i.
   * \param[in] val_windgust_j - Wind gust of the point j.
   */
  inline void SetWindGust(su2double *val_windgust_i, su2double *val_windgust_j) {
    WindGust_i = val_windgust_i;
    WindGust_j = val_windgust_j;
  }

  /*!
   * \brief Set the wind gust derivatives values.
   * \param[in] val_windgust_i - Wind gust derivatives of the point i.
   * \param[in] val_windgust_j - Wind gust derivatives of the point j.
   */
  inline void SetWindGustDer(su2double *val_windgustder_i, su2double *val_windgustder_j) {
    WindGustDer_i = val_windgustder_i;
    WindGustDer_j = val_windgustder_j;
  }

  /*!
   * \brief Set the value of the pressure.
   * \param[in] val_pressure_i - Value of the pressure at point i.
   * \param[in] val_pressure_j - Value of the pressure at point j.
   */
  void SetPressure(su2double val_pressure_i, su2double val_pressure_j) {
    Pressure_i = val_pressure_i;
    Pressure_j = val_pressure_j;
  }

  /*!
   * \brief Set the value of the density for the incompressible solver.
   * \param[in] val_densityinc_i - Value of the pressure at point i.
   * \param[in] val_densityinc_j - Value of the pressure at point j.
   */
  void SetDensity(su2double val_densityinc_i, su2double val_densityinc_j) {
    DensityInc_i = val_densityinc_i;
    DensityInc_j = val_densityinc_j;
  }

  /*!
   * \brief Set the value of the beta for incompressible flows.
   * \param[in] val_betainc2_i - Value of beta for incompressible flows at point i.
   * \param[in] val_betainc2_j - Value of beta for incompressible flows at point j.
   */
  inline void SetBetaInc2(su2double val_betainc2_i, su2double val_betainc2_j) {
    BetaInc2_i = val_betainc2_i;
    BetaInc2_j = val_betainc2_j;
  }

  /*!
   * \brief Set the value of the sound speed.
   * \param[in] val_soundspeed_i - Value of the sound speed at point i.
   * \param[in] val_soundspeed_j - Value of the sound speed at point j.
   */
  inline void SetSoundSpeed(su2double val_soundspeed_i, su2double val_soundspeed_j) {
    SoundSpeed_i = val_soundspeed_i;
    SoundSpeed_j = val_soundspeed_j;
  }

  /*!
   * \brief Set the value of the temperature.
   * \param[in] val_temp_i - Value of the temperature at point i.
   * \param[in] val_temp_j - Value of the temperature at point j.
   */
  inline void SetTemperature(su2double val_temp_i, su2double val_temp_j) {
    Temp_i = val_temp_i;
    Temp_j = val_temp_j;
  }

  /*!
   * \brief Set the value of the enthalpy.
   * \param[in] val_enthalpy_i - Value of the enthalpy at point i.
   * \param[in] val_enthalpy_j - Value of the enthalpy at point j.
   */
  inline void SetEnthalpy(su2double val_enthalpy_i, su2double val_enthalpy_j) {
    Enthalpy_i = val_enthalpy_i;
    Enthalpy_j = val_enthalpy_j;
  }

  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] val_lambda_i - Value of the spectral radius at point i.
   * \param[in] val_lambda_j - Value of the spectral radius at point j.
   */
  inline void SetLambda(su2double val_lambda_i, su2double val_lambda_j) {
    Lambda_i = val_lambda_i;
    Lambda_j = val_lambda_j;
  }

  /*!
   * \brief Set the value of undivided laplacian.
   * \param[in] val_und_lapl_i Undivided laplacian at point i.
   * \param[in] val_und_lapl_j Undivided laplacian at point j.
   */
  inline void SetUndivided_Laplacian(su2double *val_und_lapl_i, su2double *val_und_lapl_j) {
    Und_Lapl_i = val_und_lapl_i;
    Und_Lapl_j = val_und_lapl_j;
  }

  /*!
   * \brief Set the value of the pressure sensor.
   * \param[in] val_sensor_i Pressure sensor at point i.
   * \param[in] val_sensor_j Pressure sensor at point j.
   */
  inline void SetSensor(su2double val_sensor_i, su2double val_sensor_j) {
    Sensor_i = val_sensor_i;
    Sensor_j = val_sensor_j;
  }

  /*!
   * \brief Set the number of neighbor to a point.
   * \param[in] val_neighbor_i - Number of neighbor to point i.
   * \param[in] val_neighbor_j - Number of neighbor to point j.
   */
  inline void SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
    Neighbor_i = val_neighbor_i;
    Neighbor_j = val_neighbor_j;
  }

  /*!
   * \brief Set the value of the normal vector to the face between two points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  inline void SetNormal(su2double *val_normal) { Normal = val_normal; }

  /*!
   * \brief Set the value of the volume of the control volume.
   * \param[in] val_volume Volume of the control volume.
   */
  inline void SetVolume(su2double val_volume) { Volume = val_volume; }

  /*!
  * \brief Sets the values of the roe dissipation.
  * \param[in] diss_i - Dissipation value at node i
  * \param[in] diss_j - Dissipation value at node j
  */
  inline void SetDissipation(su2double diss_i, su2double diss_j) {
    Dissipation_i = diss_i;
    Dissipation_j = diss_j;
  }

  /*!
  * \brief Get the final Roe dissipation factor.
  */
  inline su2double GetDissipation() const { return Dissipation_ij; }

  /*!
   * \brief Get the inviscid fluxes.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_pressure - Value of the pressure.
   * \param[in] val_enthalpy - Value of the enthalpy.
   */
  void GetInviscidFlux(su2double val_density, su2double *val_velocity, su2double val_pressure, su2double val_enthalpy);

  /*!
   * \brief Compute the projected inviscid flux vector.
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_enthalpy - Pointer to the enthalpy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidProjFlux(su2double *val_density, su2double *val_velocity,
                           su2double *val_pressure, su2double *val_enthalpy,
                           su2double *val_normal, su2double *val_Proj_Flux);

  /*!
   * \brief Compute the projected inviscid flux vector for incompresible simulations
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidIncProjFlux(su2double *val_density, su2double *val_velocity,
                                  su2double *val_pressure, su2double *val_betainc2,
                                  su2double *val_enthalpy,
                                  su2double *val_normal, su2double *val_Proj_Flux);

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices.
   * \param[in] val_velocity Pointer to the velocity.
   * \param[in] val_energy Value of the energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(su2double *val_velocity, su2double *val_energy,
                          su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_tensor);

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices (incompressible).
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidIncProjJac(su2double *val_density, su2double *val_velocity,
                             su2double *val_betainc2, su2double *val_normal,
                             su2double val_scale,
                             su2double **val_Proj_Jac_tensor);

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices (overload for low speed preconditioner version).
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_cp - Value of the specific heat at constant pressure.
   * \param[in] val_temperature - Value of the temperature.
   * \param[in] val_dRhodT - Value of the derivative of density w.r.t. temperature.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidIncProjJac(su2double *val_density,
                             su2double *val_velocity,
                             su2double *val_betainc2,
                             su2double *val_cp,
                             su2double *val_temperature,
                             su2double *val_dRhodT,
                             su2double *val_normal,
                             su2double val_scale,
                             su2double **val_Proj_Jac_Tensor);

  /*!
   * \brief Compute the low speed preconditioning matrix.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_cp - Value of the specific heat at constant pressure.
   * \param[in] val_temperature - Value of the temperature.
   * \param[in] val_dRhodT - Value of the derivative of density w.r.t. temperature.
   * \param[out] val_Precon - Pointer to the preconditioning matrix.
   */
  void GetPreconditioner(su2double *val_density,
                         su2double *val_velocity,
                         su2double *val_betainc2,
                         su2double *val_cp,
                         su2double *val_temperature,
                         su2double *val_drhodt,
                         su2double **val_Precon);

  /*!
   * \brief Compute the projection of the preconditioned inviscid Jacobian matrices.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetPreconditionedProjJac(su2double *val_density,
                                su2double *val_velocity,
                                su2double *val_betainc2,
                                su2double *val_normal,
                                su2double **val_Proj_Jac_Tensor);

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices for general fluid model.
   * \param[in] val_velocity Pointer to the velocity.
   * \param[in] val_energy Value of the energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(su2double *val_velocity, su2double *val_enthalphy,
                          su2double *val_chi, su2double *val_kappa,
                          su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_tensor);

  /*!
   * \brief Mapping between primitives variables P and conservatives variables C.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_PrimVar - Mean Value of the secondary variables.
   * \param[out] val_Jac_PC - Pointer to the Jacobian dPdC.
   */
  void GetPrimitive2Conservative (su2double *val_Mean_PrimVar,
                                  su2double *val_Mean_SecVar,
                                  su2double **val_Jac_PC);

  /*!
   * \overload
   * \brief Computation of the matrix P for a generic fluid model
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_enthalpy - Value of the Enthalpy
   * \param[in] val_chi - Value of the derivative of Pressure with respect to the Density.
   * \param[in] val_kappa - Value of the derivative of Pressure with respect to the volume specific Static Energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPMatrix(su2double *val_density, su2double *val_velocity,
                  su2double *val_soundspeed, su2double *val_enthalpy,
                  su2double *val_chi, su2double *val_kappa,
                  su2double *val_normal, su2double **val_p_tensor);

  /*!
   * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in
   *        the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPMatrix(su2double *val_density, su2double *val_velocity,
                  su2double *val_soundspeed, su2double *val_normal,
                  su2double **val_p_tensor);

  /*!
   * \brief Computation of the matrix Rinv*Pe.
   * \param[in] Beta2 - A variable in used to define Pe matrix.
   * \param[in] val_enthalpy - value of the enthalpy.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[in] val_velocity - value of the velocity.
   * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetinvRinvPe(su2double Beta2, su2double val_enthalpy, su2double val_soundspeed,
                    su2double val_density, su2double* val_velocity,
                    su2double** val_invR_invPe);

  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_pressure - value of the pressure.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[in] val_velocity - value of the velocity.
   * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetRMatrix(su2double val_pressure, su2double val_soundspeed,
                  su2double val_density, su2double* val_velocity,
                  su2double** val_invR_invPe);
  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[out] R_Matrix - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetRMatrix(su2double val_soundspeed, su2double val_density, su2double **R_Matrix);

  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[out] L_Matrix - Pointer to the matrix of conversion from conserved to entropic variables.
   */
  void GetLMatrix(su2double val_soundspeed, su2double val_density, su2double **L_Matrix);

  /*!
   * \brief Computation of the flow Residual Jacoboan Matrix for Non Reflecting BC.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[out] R_c - Residual Jacoboan Matrix
   * \param[out] R_c_inv- inverse of the Residual Jacoboan Matrix .
   */
  void ComputeResJacobianGiles(CFluidModel *FluidModel, su2double pressure, su2double density, su2double *turboVel,
                               su2double alphaInBC, su2double gammaInBC,  su2double **R_c, su2double **R_c_inv);

  /*!
   * \brief Computate the inverse of a 3x3 matrix
   * \param[in]  matrix     - the matrix to invert
   * \param[out] invMatrix  - inverse matrix.
   */
  void InvMatrix3D(su2double **matrix, su2double **invMatrix);

  /*!
   * \brief Computate the inverse of a 4x4 matrix
   * \param[in]  matrix    - the matrix to invert
   * \param[out] invMatrix - inverse matrix.
   */
  void InvMatrix4D(su2double **matrix, su2double **invMatrix);

  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[in] prim_jump - pointer to the vector containing the primitive variable jump (drho, dV, dp).
   * \param[out] char_jump - pointer to the vector containing the characteristic variable jump.
   */
  void GetCharJump(su2double val_soundspeed, su2double val_density, su2double *prim_jump, su2double *char_jump);

  /*!
   * \brief Computation of the matrix Td, this matrix diagonalize the preconditioned conservative Jacobians
   *        in the form $Tg |Lambda| Td = Pc{-1}|Pc (A.Normal)|$.
   * \param[in] Beta2 - A variable in used to define absPeJacobian matrix.
   * \param[in] r_hat - A variable in used to define absPeJacobian matrix.
   * \param[in] s_hat - A variable in used to define absPeJacobian matrix.
   * \param[in] t_hat - A variable in used to define absPeJacobian matrix.
   * \param[in] rB2a2 - A variable in used to define absPeJacobian matrix.
   * \param[in] val_Lambda - Eigenvalues of the Preconditioned Jacobian.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_absPeJac - Pointer to the Preconditioned Jacobian matrix.
   */
  void GetPrecondJacobian(su2double Beta2, su2double r_hat, su2double s_hat, su2double t_hat,
                          su2double rB2a2, su2double* val_Lambda, su2double* val_normal, su2double** val_absPeJac);

  /*!
   * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
   * in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(su2double **val_invp_tensor, su2double *val_density,
                      su2double *val_velocity, su2double *val_soundspeed,
                      su2double *val_chi, su2double *val_kappa,
                      su2double *val_normal);

  /*!
   * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(su2double *val_density, su2double *val_velocity,
                      su2double *val_soundspeed, su2double *val_normal,
                      su2double **val_invp_tensor);

  /*!
   * \brief Compute viscous residual and jacobian.
   */
  void GetAdjViscousFlux_Jac(su2double Pressure_i, su2double Pressure_j, su2double Density_i, su2double Density_j,
                             su2double ViscDens_i, su2double ViscDens_j, su2double *Velocity_i, su2double *Velocity_j,
                             su2double sq_vel_i, su2double sq_vel_j,
                             su2double XiDens_i, su2double XiDens_j, su2double **Mean_GradPhi, su2double *Mean_GradPsiE,
                             su2double dPhiE_dn, su2double *Normal, su2double *Edge_Vector, su2double dist_ij_2, su2double *val_residual_i,
                             su2double *val_residual_j,
                             su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji,
                             su2double **val_Jacobian_jj, bool implicit);

  /*!
   * \brief Computation of the projected inviscid lambda (eingenvalues).
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_Lambda_Vector - Pointer to Lambda matrix.
   */
  void GetJacInviscidLambda_fabs(su2double *val_velocity, su2double val_soundspeed,
                                 su2double *val_normal, su2double *val_Lambda_Vector);

  /*!
   * \brief Compute the numerical residual.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ComputeResidual(su2double *val_residual, CConfig* config) { }

  /*!
   * \overload
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                      su2double **val_Jacobian_j, CConfig* config) { }

  /*!
   * \overload For numerics classes that store the residual/flux and Jacobians internally.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  inline virtual ResidualType<> ComputeResidual(const CConfig* config) { return ResidualType<>(nullptr,nullptr,nullptr); }

  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ComputeResidual(su2double *val_residual_i,
                                      su2double *val_residual_j,
                                      su2double **val_Jacobian_ii,
                                      su2double **val_Jacobian_ij,
                                      su2double **val_Jacobian_ji,
                                      su2double **val_Jacobian_jj, CConfig* config) { }

  /*!
   * \overload
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ComputeResidual(su2double *val_resconv_i, su2double *val_resvisc_i,
                                      su2double *val_resconv_j, su2double *val_resvisc_j,
                                      su2double **val_Jacobian_ii,
                                      su2double **val_Jacobian_ij,
                                      su2double **val_Jacobian_ji,
                                      su2double **val_Jacobian_jj, CConfig* config) { }

  /*!
   * \overload
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  inline virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig* config) { }

  /*!
   * \brief Residual for transition problems.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   * \param[out] gamma_sep
   */
  inline virtual void ComputeResidual_TransLM(su2double *val_residual,
                                              su2double **val_Jacobian_i,
                                              su2double **val_Jacobian_j, CConfig* config,
                                              su2double &gamma_sep) { }

  /*!
   * \brief Set intermittency for numerics (used in SA with LM transition model)
   */
  inline virtual void SetIntermittency(su2double intermittency_in) { }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  inline virtual void SetProduction(su2double val_production) { }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  inline virtual void SetDestruction(su2double val_destruction) { }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  inline virtual void SetCrossProduction(su2double val_crossproduction) { }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  inline virtual su2double GetProduction(void) const { return 0; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  inline virtual su2double GetDestruction(void) const { return 0; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  inline virtual su2double GetCrossProduction(void) const { return 0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetGammaBC(void) const { return 0.0; }

  /*!
   * \brief A virtual member to compute the tangent matrix in structural problems
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  inline virtual void Compute_Tangent_Matrix(CElement *element_container, const CConfig* config) { }

  /*!
   * \brief A virtual member to compute the nodal stress term in non-linear structural problems
   * \param[in] element_container - Definition of the particular element integrated.
   */
  inline virtual void Compute_NodalStress_Term(CElement *element_container, const CConfig* config) { }

  /*!
   * \brief Set the element-based local Young's modulus in mesh problems
   * \param[in] iElem - Element index.
   * \param[in] val_E - Value of elasticity modulus.
   */
  inline virtual void SetMeshElasticProperties(unsigned long iElem, su2double val_E) { }

  /*!
   * \brief A virtual member to set the value of the design variables
   * \param[in] i_DV - Index of the design variable.
   * \param[in] val_DV - Value of the design variable
   */
  inline virtual void Set_DV_Val(unsigned short i_DV, su2double val_DV) { }

  /*!
   * \brief A virtual member to retrieve the value of the design variables
   * \param[in] i_DV - Index of the design variable.
   */
  inline virtual su2double Get_DV_Val(unsigned short i_DV) const { return 0.0; }

  /*!
   * \brief A virtual member to set the electric field
   * \param[in] EField_DV - New electric field computed by adjoint methods.
   */
  inline virtual void Set_ElectricField(unsigned short i_DV, su2double val_EField) { }

  /*!
   * \brief A virtual member to set the material properties
   * \param[in] iVal - Index of the region of concern
   * \param[in] val_E - Value of the Young Modulus.
   * \param[in] val_Nu - Value of the Poisson's ratio.
   */
  inline virtual void SetMaterial_Properties(unsigned short iVal, su2double val_E, su2double val_Nu) { }

  /*!
   * \brief A virtual member to set the material properties
   * \param[in] iVal - Index of the region of concern
   * \param[in] val_Rho - Value of the density (inertial effects).
   * \param[in] val_Rho_DL - Value of the density (dead load effects).
   */
  inline virtual void SetMaterial_Density(unsigned short iVal, su2double val_Rho, su2double val_Rho_DL) { }

  /*!
   * \brief A virtual member to compute the mass matrix
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  inline virtual void Compute_Mass_Matrix(CElement *element_container, const CConfig* config) { }

  /*!
   * \brief A virtual member to compute the residual component due to dead loads
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  inline virtual void Compute_Dead_Load(CElement *element_container, const CConfig* config) { }

  /*!
   * \brief A virtual member to compute the averaged nodal stresses
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  inline virtual void Compute_Averaged_NodalStress(CElement *element_container, const CConfig* config) { }

  /*!
   * \brief Computes a basis of orthogonal vectors from a supplied vector
   * \param[in] config - Normal vector
   */
  void CreateBasis(su2double *val_Normal);

  /*!
   * \brief Set the value of the Tauwall
   * \param[in] val_tauwall_i - Tauwall at point i
   * \param[in] val_tauwall_j - Tauwall at point j
   */
  inline virtual void SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) { }

  /*!
   * \brief Set the value of the bollean flag to use (or not) the wall shear stress from the wall function.
   * \param[in] val_tauwallflag_i - Flag for Tauwall at point i
   * \param[in] val_tauwallflag_j - Flag for Tauwall at point j
   */
  inline virtual void SetTauWall_Flag(bool val_tauwallflag_i, bool val_tauwallflag_j) { }

  /*!
   * \brief - Calculate the central/upwind blending function for a face
   *
   * At its simplest level, this function will calculate the average
   * value of the "Roe dissipation" or central/upwind blending function
   * at a face. If Ducros shock sensors are used, then this method will
   * also adjust the central/upwind blending to account for shocks.
   *
   * \param[in] Dissipation_i - Value of the blending function at point i
   * \param[in] Dissipation_j - Value of the bledning function at point j
   * \param[in] Sensor_i - Shock sensor at point i
   * \param[in] Sensor_j - Shock sensor at point j
   * \param[in] config - The configuration of the problem
   * \return Dissipation_ij - Blending parameter at face
   */
  su2double GetRoe_Dissipation(const su2double Dissipation_i, const su2double Dissipation_j,
                               const su2double Sensor_i, const su2double Sensor_j, const CConfig* config) const;

  /*!
   * \brief Set the value of the radiation variable.
   * \param[in] val_radvar_i - Value of the turbulent variable at point i.
   * \param[in] val_radvar_j - Value of the turbulent variable at point j.
   */
  inline virtual void SetRadVar(const su2double *val_radvar_i, const su2double *val_radvar_j) { }

  /*!
   * \brief Set the gradient of the radiation variables.
   * \param[in] val_radvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_radvar_grad_j - Gradient of the turbulent variable at point j.
   */
  inline virtual void SetRadVarGradient(const su2double* const* val_radvar_grad_i, const su2double* const* val_radvar_grad_j) { }

  /*!
   * \brief Set the gradient of the radiation variables.
   * \param[in] val_radvar_source - Source term (and jacobian term) of the radiative heat transfer.
   */
  inline void SetRadVarSource(const su2double *val_radvar_source) { RadVar_Source = val_radvar_source; }

  /*!
   * \brief Decomposes the symmetric matrix A_ij, into eigenvectors and eigenvalues
   * \param[in] A_i: symmetric matrix to be decomposed
   * \param[in] Eig_Vec: strores the eigenvectors
   * \param[in] Eig_Val: stores the eigenvalues
   * \param[in] n: order of matrix A_ij
   */
  static void EigenDecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val, unsigned short n);

  /*!
   * \brief Recomposes the eigenvectors and eigenvalues into a matrix
   * \param[in] A_ij: recomposed matrix
   * \param[in] Eig_Vec: eigenvectors
   * \param[in] Eig_Val: eigenvalues
   * \param[in] n: order of matrix A_ij
   */
  static void EigenRecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val, unsigned short n);

  /*!
   * \brief tred2
   * \param[in] V: matrix that needs to be decomposed
   * \param[in] d: holds eigenvalues
   * \param[in] e: supplemental data structure
   * \param[in] n: order of matrix V
   */
  static void tred2(su2double **V, su2double *d, su2double *e, unsigned short n);

  /*!
   * \brief tql2
   * \param[in] V: matrix that will hold the eigenvectors
   * \param[in] d: array that will hold the ordered eigenvalues
   * \param[in] e: supplemental data structure
   * \param[in] n: order of matrix V
   */
  static void tql2(su2double **V, su2double *d, su2double *e, unsigned short n);

};

/*!
 * /brief Alias for a "do nothing" source term class
 */
using CSourceNothing = CNumerics;
