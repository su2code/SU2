/*!
 * \file CNumerics.hpp
 * \brief Declaration of the base numerics class, the
 *        implementation is in the CNumerics.cpp file.
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

#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"

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
  su2double Prandtl_Lam;      /*!< \brief Laminar Prandtl's number. */
  su2double Prandtl_Turb;     /*!< \brief Turbulent Prandtl's number. */
  su2double MassFlux;         /*!< \brief Mass flux across edge. */
  su2double
  *Proj_Flux_Tensor;  /*!< \brief Flux tensor projected in a direction. */
  su2double **tau;    /*!< \brief Viscous stress tensor. */
  const su2double delta [3][3] = {{1.0, 0.0, 0.0},{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}; /*!< \brief Identity matrix. */
  const su2double
  *Diffusion_Coeff_i, /*!< \brief Species diffusion coefficients at point i. */
  *Diffusion_Coeff_j; /*!< \brief Species diffusion coefficients at point j. */
  su2double
  Laminar_Viscosity_i,   /*!< \brief Laminar viscosity at point i. */
  Laminar_Viscosity_j;   /*!< \brief Laminar viscosity at point j. */
  su2double
  Thermal_Conductivity_i,    /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_j,    /*!< \brief Thermal conductivity at point j. */
  Thermal_Conductivity_ve_i, /*!< \brief vibrational-electronic Thermal conductivity at point i. */
  Thermal_Conductivity_ve_j, /*!< \brief vibrational-electronic Thermal conductivity at point j. */
  Cp_i,               /*!< \brief Cp at point i. */
  Cp_j;               /*!< \brief Cp at point j. */
  su2double
  Eddy_Viscosity_i,  /*!< \brief Eddy viscosity at point i. */
  Eddy_Viscosity_j;  /*!< \brief Eddy viscosity at point j. */
  su2double
  turb_ke_i,  /*!< \brief Turbulent kinetic energy at point i. */
  turb_ke_j;  /*!< \brief Turbulent kinetic energy at point j. */
  su2double
  intermittency_eff_i, /*!< \brief effective intermittency at point i. */
  intermittency_i; /*!< \brief intermittency at point i. */
  su2double
  Pressure_i,  /*!< \brief Pressure at point i. */
  Pressure_j;  /*!< \brief Pressure at point j. */
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
  SoundSpeed_i,  /*!< \brief Sound speed at point i. */
  SoundSpeed_j;  /*!< \brief Sound speed at point j. */
  su2double
  Enthalpy_i,  /*!< \brief Enthalpy at point i. */
  Enthalpy_j;  /*!< \brief Enthalpy at point j. */
  su2double
  dist_i,  /*!< \brief Distance of point i to the nearest wall. */
  dist_j;  /*!< \brief Distance of point j to the nearest wall. */
  const su2double
  *Und_Lapl_i,  /*!< \brief Undivided laplacians at point i. */
  *Und_Lapl_j;  /*!< \brief Undivided laplacians at point j. */
  su2double
  Sensor_i,  /*!< \brief Pressure sensor at point i. */
  Sensor_j;  /*!< \brief Pressure sensor at point j. */
  const su2double
  *GridVel_i,  /*!< \brief Grid velocity at point i. */
  *GridVel_j;  /*!< \brief Grid velocity at point j. */
  const su2double
  *U_i,           /*!< \brief Vector of conservative variables at point i. */
  *U_j;           /*!< \brief Vector of conservative variables at point j. */
  const su2double
  *V_i,     /*!< \brief Vector of primitive variables at point i. */
  *V_j;     /*!< \brief Vector of primitive variables at point j. */
  const su2double
  *S_i,     /*!< \brief Vector of secondary variables at point i. */
  *S_j;     /*!< \brief Vector of secondary variables at point j. */
  const su2double
  *Psi_i,    /*!< \brief Vector of adjoint variables at point i. */
  *Psi_j;    /*!< \brief Vector of adjoint variables at point j. */
  const su2double
  *ScalarVar_i,   /*!< \brief Vector of scalar variables at point i. */
  *ScalarVar_j;   /*!< \brief Vector of scalar variables at point j. */
  const su2double
  *TransVar_i,  /*!< \brief Vector of turbulent variables at point i. */
  *TransVar_j;  /*!< \brief Vector of turbulent variables at point j. */
  const su2double
  *TurbPsi_i,  /*!< \brief Vector of adjoint turbulent variables at point i. */
  *TurbPsi_j;  /*!< \brief Vector of adjoint turbulent variables at point j. */
  CMatrixView<const su2double>
  ConsVar_Grad_i,  /*!< \brief Gradient of conservative variables at point i. */
  ConsVar_Grad_j,  /*!< \brief Gradient of conservative variables at point j. */
  ConsVar_Grad,    /*!< \brief Gradient of conservative variables which is a scalar. */
  PrimVar_Grad_i,  /*!< \brief Gradient of primitive variables at point i. */
  PrimVar_Grad_j,  /*!< \brief Gradient of primitive variables at point j. */
  PsiVar_Grad_i,   /*!< \brief Gradient of adjoint variables at point i. */
  PsiVar_Grad_j,   /*!< \brief Gradient of adjoint variables at point j. */
  ScalarVar_Grad_i,  /*!< \brief Gradient of scalar variables at point i. */
  ScalarVar_Grad_j,  /*!< \brief Gradient of scalar variables at point j. */
  TransVar_Grad_i, /*!< \brief Gradient of turbulent variables at point i. */
  TransVar_Grad_j, /*!< \brief Gradient of turbulent variables at point j. */
  TurbPsi_Grad_i,  /*!< \brief Gradient of adjoint turbulent variables at point i. */
  TurbPsi_Grad_j,  /*!< \brief Gradient of adjoint turbulent variables at point j. */
  AuxVar_Grad_i,   /*!< \brief Gradient of an auxiliary variable at point i. */
  AuxVar_Grad_j;   /*!< \brief Gradient of an auxiliary variable at point i. */
  su2double
  LocalGridLength_i; /*!< \brief Local grid length at point i. */
  const su2double *RadVar_Source;  /*!< \brief Source term from the radiative heat transfer equation. */
  const su2double
  *Coord_i,      /*!< \brief Cartesians coordinates of point i. */
  *Coord_j;      /*!< \brief Cartesians coordinates of point j. */
  unsigned short
  Neighbor_i,  /*!< \brief Number of neighbors of the point i. */
  Neighbor_j;  /*!< \brief Number of neighbors of the point j. */
  const su2double *Normal = nullptr;      /*!< \brief Normal vector, its norm is the area of the face. */
  su2double UnitNormal[MAXNDIM] = {0.0};  /*!< \brief Unitary normal vector. */
  su2double
  TimeStep,    /*!< \brief Time step useful in dual time method. */
  Area,        /*!< \brief Area of the face i-j. */
  Volume,      /*!< \brief Volume of the control volume around point i. */
  AvgVolume;    /*!< \brief Average of the control Volume around point i for vorticity confinement parameter correction */
  su2double vel2_inf;     /*!< \brief value of the square of freestream speed. */
  const su2double
  *WindGust_i,  /*!< \brief Wind gust at point i. */
  *WindGust_j;  /*!< \brief Wind gust at point j. */
  const su2double *Vorticity_i, *Vorticity_j;    /*!< \brief Vorticity. */
  su2double StrainMag_i, StrainMag_j;      /*!< \brief Strain rate magnitude. */
  su2double Dissipation_i, Dissipation_j;  /*!< \brief Dissipation. */
  su2double Dissipation_ij;
  su2double roughness_i = 0.0,             /*!< \brief Roughness of the wall nearest to point i. */
  roughness_j = 0.0;                       /*!< \brief Roughness of the wall nearest to point j. */

  su2double MeanPerturbedRSM[3][3];   /*!< \brief Perturbed Reynolds stress tensor  */
  SST_ParsedOptions sstParsedOptions; /*!< \brief additional options for the SST turbulence model */
  unsigned short Eig_Val_Comp;    /*!< \brief Component towards which perturbation is perfromed */
  su2double uq_delta_b;           /*!< \brief Magnitude of perturbation */
  su2double uq_urlx;              /*!< \brief Under-relaxation factor for numerical stability */
  bool uq_permute;                /*!< \brief Flag for eigenvector permutation */

  bool nemo;                      /*!< \brief Flag for NEMO problems  */

  bool bounded_scalar = false;    /*!< \brief Flag for bounded scalar problem */

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

    su2double operator[] (unsigned long idx) const { return residual[idx]; }
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
  void SetVorticity(const su2double *val_vorticity_i, const su2double *val_vorticity_j) {
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
  inline void SetConservative(const su2double *val_u_i, const su2double *val_u_j) {
    U_i = val_u_i;
    U_j = val_u_j;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  inline void SetPrimitive(const su2double *val_v_i, const su2double *val_v_j) {
    V_i = val_v_i;
    V_j = val_v_j;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  inline void SetSecondary(const su2double *val_s_i, const su2double *val_s_j) {
    S_i = val_s_i;
    S_j = val_s_j;
  }

  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_i - Gradient of the conservative variable at point i.
   * \param[in] val_consvar_grad_j - Gradient of the conservative variable at point j.
   */
  inline void SetConsVarGradient(CMatrixView<const su2double> val_consvar_grad_i,
                                 CMatrixView<const su2double> val_consvar_grad_j) {
    ConsVar_Grad_i = val_consvar_grad_i;
    ConsVar_Grad_j = val_consvar_grad_j;
  }

  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad - Gradient of the conservative variable which is a scalar.
   */
  inline void SetConsVarGradient(CMatrixView<const su2double> val_consvar_grad) { ConsVar_Grad = val_consvar_grad; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
   * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
   */
  void SetPrimVarGradient(CMatrixView<const su2double> val_primvar_grad_i,
                          CMatrixView<const su2double> val_primvar_grad_j) {
    PrimVar_Grad_i = val_primvar_grad_i;
    PrimVar_Grad_j = val_primvar_grad_j;
  }

  /*!
   * \brief Set the value of the adjoint variable.
   * \param[in] val_psi_i - Value of the adjoint variable at point i.
   * \param[in] val_psi_j - Value of the adjoint variable at point j.
   */
  inline void SetAdjointVar(const su2double *val_psi_i, const su2double *val_psi_j) {
    Psi_i = val_psi_i;
    Psi_j = val_psi_j;
  }

  /*!
   * \brief Set the gradient of the adjoint variables.
   * \param[in] val_psivar_grad_i - Gradient of the adjoint variable at point i.
   * \param[in] val_psivar_grad_j - Gradient of the adjoint variable at point j.
   */
  inline void SetAdjointVarGradient(CMatrixView<const su2double> val_psivar_grad_i,
                                    CMatrixView<const su2double> val_psivar_grad_j) {
    PsiVar_Grad_i = val_psivar_grad_i;
    PsiVar_Grad_j = val_psivar_grad_j;
  }

  /*!
   * \brief Set the value of the scalar variable.
   * \param[in] val_scalarvar_i - Value of the scalar variable at point i.
   * \param[in] val_scalarvar_j - Value of the scalar variable at point j.
   */
  inline void SetScalarVar(const su2double *val_scalarvar_i, const su2double *val_scalarvar_j) {
    ScalarVar_i = val_scalarvar_i;
    ScalarVar_j = val_scalarvar_j;
  }

  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_transvar_i - Value of the turbulent variable at point i.
   * \param[in] val_transvar_j - Value of the turbulent variable at point j.
   */
  inline void SetTransVar(const su2double *val_transvar_i, const su2double *val_transvar_j) {
    TransVar_i = val_transvar_i;
    TransVar_j = val_transvar_j;
  }

  /*!
   * \brief Set the gradient of the scalar variables.
   * \param[in] val_scalarvar_grad_i - Gradient of the scalar variable at point i.
   * \param[in] val_scalarvar_grad_j - Gradient of the scalar variable at point j.
   */
  inline void SetScalarVarGradient(CMatrixView<const su2double> val_scalarvar_grad_i,
                                   CMatrixView<const su2double> val_scalarvar_grad_j) {
    ScalarVar_Grad_i = val_scalarvar_grad_i;
    ScalarVar_Grad_j = val_scalarvar_grad_j;
  }

  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  inline void SetTransVarGradient(CMatrixView<const su2double> val_transvar_grad_i,
                                  CMatrixView<const su2double> val_transvar_grad_j) {
    TransVar_Grad_i = val_transvar_grad_i;
    TransVar_Grad_j = val_transvar_grad_j;
  }

  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_transvar_i - Value of the turbulent variable at point i.
   */
  inline void SetLocalGridLength(const su2double val_localGridLength_i) {
    LocalGridLength_i = val_localGridLength_i;
  }

  /*!
   * \brief Set the value of the adjoint turbulent variable.
   * \param[in] val_turbpsivar_i - Value of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_j - Value of the adjoint turbulent variable at point j.
   */
  inline void SetTurbAdjointVar(const su2double *val_turbpsivar_i, const su2double *val_turbpsivar_j) {
    TurbPsi_i = val_turbpsivar_i;
    TurbPsi_j = val_turbpsivar_j;
  }

  /*!
   * \brief Set the gradient of the adjoint turbulent variables.
   * \param[in] val_turbpsivar_grad_i - Gradient of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_grad_j - Gradient of the adjoint turbulent variable at point j.
   */
  inline void SetTurbAdjointGradient(CMatrixView<const su2double> val_turbpsivar_grad_i,
                                     CMatrixView<const su2double> val_turbpsivar_grad_j) {
    TurbPsi_Grad_i = val_turbpsivar_grad_i;
    TurbPsi_Grad_j = val_turbpsivar_grad_j;
  }

  /*!
   * \brief Compute the mean rate of strain matrix.
   * \details The parameter primvargrad can be e.g. PrimVar_Grad_i or Mean_GradPrimVar.
   * \param[in] nDim - 2 or 3
   * \param[out] rateofstrain - Rate of strain matrix
   * \param[in] velgrad - A velocity gradient matrix.
   * \tparam Mat1 - any type that supports the [][] interface
   * \tparam Mat2 - any type that supports the [][] interface
   */
  template<class Mat1, class Mat2>
  FORCEINLINE static void ComputeMeanRateOfStrainMatrix(size_t nDim, Mat1& rateofstrain, const Mat2& velgrad){

    /* --- Calculate the rate of strain tensor, using mean velocity gradients --- */

    if (nDim == 3) {
      rateofstrain[0][0] = velgrad[0][0];
      rateofstrain[1][1] = velgrad[1][1];
      rateofstrain[2][2] = velgrad[2][2];
      rateofstrain[0][1] = 0.5 * (velgrad[0][1] + velgrad[1][0]);
      rateofstrain[0][2] = 0.5 * (velgrad[0][2] + velgrad[2][0]);
      rateofstrain[1][2] = 0.5 * (velgrad[1][2] + velgrad[2][1]);
      rateofstrain[1][0] = rateofstrain[0][1];
      rateofstrain[2][1] = rateofstrain[1][2];
      rateofstrain[2][0] = rateofstrain[0][2];
    }
    else { // nDim==2
      rateofstrain[0][0] = velgrad[0][0];
      rateofstrain[1][1] = velgrad[1][1];
      rateofstrain[2][2] = 0.0;
      rateofstrain[0][1] = 0.5 * (velgrad[0][1] + velgrad[1][0]);
      rateofstrain[0][2] = 0.0;
      rateofstrain[1][2] = 0.0;
      rateofstrain[1][0] = rateofstrain[0][1];
      rateofstrain[2][1] = rateofstrain[1][2];
      rateofstrain[2][0] = rateofstrain[0][2];
    }
  }

  /*!
   * \brief Compute the stress tensor from the velocity gradients.
   * \details To obtain the Reynolds stress tensor +(u_i' u_j')~, divide the result
   * of this function by (-rho). The argument density is only used if turb_ke is not 0.
   * To select the velocity gradient components from a primitive variable gradient PrimVar_Grad_i,
   * write PrimVar_Grad_i+1.
   * If <code>nDim==2</code>, we use the same formula but only only access the entries [0][0]..[1][1] of
   * stress and velgrad. If <code>reynolds3x3</code> is true, the other non-diagonal entries of stress
   * set to zero, and <code>stress[2][2]</code> to some value.
   * \param[in] nDim - Dimension of the flow problem, 2 or 3
   * \param[out] stress - Stress tensor
   * \param[in] velgrad - A velocity gradient matrix.
   * \param[in] viscosity - Viscosity
   * \param[in] density - Density
   * \param[in] turb_ke - Turbulent kinetic energy, for the turbulent stress tensor
   * \param[in] reynolds3x3 - If true, write to the third row and column of stress even if nDim==2.
   * \tparam Mat1 - any type that supports the [][] interface
   * \tparam Mat2 - any type that supports the [][] interface
   */
  template<class Mat1, class Mat2, class Scalar>
  FORCEINLINE static void ComputeStressTensor(size_t nDim, Mat1& stress, const Mat2& velgrad,
                                              Scalar viscosity, Scalar density=0.0,
                                              Scalar turb_ke=0.0, bool reynolds3x3=false){
    Scalar divVel = 0.0;
    for (size_t iDim = 0; iDim < nDim; iDim++) {
      divVel += velgrad[iDim][iDim];
    }
    Scalar pTerm = 2./3. * (divVel * viscosity + density * turb_ke);

    for (size_t iDim = 0; iDim < nDim; iDim++){
      for (size_t jDim = 0; jDim < nDim; jDim++){
        stress[iDim][jDim] = viscosity * (velgrad[iDim][jDim]+velgrad[jDim][iDim]);
      }
      stress[iDim][iDim] -= pTerm;
    }

    if(reynolds3x3 && nDim==2) { // fill the third row and column of Reynolds stress matrix
      stress[0][2] = stress[1][2] = stress[2][0] = stress[2][1] = 0.0;
      stress[2][2] = -pTerm;
    }
  }

  /*!
   * \brief Add a correction using a Quadratic Constitutive Relation to the stress tensor.
   *
   * See: Spalart, P. R., "Strategies for Turbulence Modelling and Simulation",
   * International Journal of Heat and Fluid Flow, Vol. 21, 2000, pp. 252-263
   *
   * \param[in] nDim: 2D or 3D.
   * \param[in] gradvel: Velocity gradients.
   * \param[in,out] tau: Shear stress tensor.
   */
  template <class Mat1, class Mat2>
  FORCEINLINE static void AddQCR(size_t nDim, const Mat1& gradvel, Mat2& tau) {
    using Scalar = typename std::decay<decltype(gradvel[0][0])>::type;

    const Scalar c_cr1 = 0.3;

    /*--- Denominator Antisymmetric normalized rotation tensor ---*/

    Scalar factor = 0.0;
    for (size_t iDim = 0; iDim < nDim; iDim++)
      for (size_t jDim = 0; jDim < nDim; jDim++)
        factor += gradvel[iDim][jDim] * gradvel[iDim][jDim];
    factor = 1.0 / sqrt(max(factor,1E-10));

    /*--- Adding the QCR contribution ---*/

    Scalar tauQCR[MAXNDIM][MAXNDIM] = {{0.0}};

    for (size_t iDim = 0; iDim < nDim; iDim++){
      for (size_t jDim = 0; jDim < nDim; jDim++){
        for (size_t kDim = 0; kDim < nDim; kDim++){
          Scalar O_ik = (gradvel[iDim][kDim] - gradvel[kDim][iDim]) * factor;
          Scalar O_jk = (gradvel[jDim][kDim] - gradvel[kDim][jDim]) * factor;
          tauQCR[iDim][jDim] += O_ik * tau[jDim][kDim] + O_jk * tau[iDim][kDim];
        }
      }
    }

    for (size_t iDim = 0; iDim < nDim; iDim++)
      for (size_t jDim = 0; jDim < nDim; jDim++)
        tau[iDim][jDim] -= c_cr1 * tauQCR[iDim][jDim];
  }

  /*!
   * \brief Perturb the Reynolds stress tensor based on parameters.
   * \param[in] nDim - Dimension of the flow problem, 2 or 3.
   * \param[in] uq_eigval_comp - Component 1C 2C 3C.
   * \param[in] uq_permute - Whether to swap order of eigen vectors.
   * \param[in] uq_delta_b - Delta_b parameter.
   * \param[in] uq_urlx - Relaxation factor.
   * \param[in] velgrad - A velocity gradient matrix.
   * \param[in] density - Density.
   * \param[in] viscosity - Eddy viscosity.
   * \param[in] turb_ke: Turbulent kinetic energy.
   * \param[out] MeanPerturbedRSM - Perturbed stress tensor.
   */
  template<class Mat1, class Mat2, class Scalar>
  NEVERINLINE static void ComputePerturbedRSM(size_t nDim, size_t uq_eigval_comp, bool uq_permute, su2double uq_delta_b,
                                              su2double uq_urlx, const Mat1& velgrad, Scalar density,
                                              Scalar viscosity, Scalar turb_ke, Mat2& MeanPerturbedRSM) {
    Scalar MeanReynoldsStress[3][3];
    ComputeStressTensor(nDim, MeanReynoldsStress, velgrad, viscosity, density, turb_ke, true);
    for (size_t iDim = 0; iDim < 3; iDim++)
      for (size_t jDim = 0; jDim < 3; jDim++)
        MeanReynoldsStress[iDim][jDim] /= -density;

    /* --- Calculate anisotropic part of Reynolds Stress tensor --- */

    Scalar A_ij[3][3];
    for (size_t iDim = 0; iDim < 3; iDim++) {
      for (size_t jDim = 0; jDim < 3; jDim++) {
        A_ij[iDim][jDim] = .5 * MeanReynoldsStress[iDim][jDim] / turb_ke;
      }
      A_ij[iDim][iDim] -= 1.0/3.0;
    }

    /* --- Get ordered eigenvectors and eigenvalues of A_ij --- */

    Scalar work[3], Eig_Vec[3][3], Eig_Val[3];
    CBlasStructure::EigenDecomposition(A_ij, Eig_Vec, Eig_Val, 3, work);

    /* compute convex combination coefficients */
    Scalar c1c = Eig_Val[2] - Eig_Val[1];
    Scalar c2c = 2.0 * (Eig_Val[1] - Eig_Val[0]);
    Scalar c3c = 3.0 * Eig_Val[0] + 1.0;

    /* define barycentric traingle corner points */
    Scalar Corners[3][2];
    Corners[0][0] = 1.0;
    Corners[0][1] = 0.0;
    Corners[1][0] = 0.0;
    Corners[1][1] = 0.0;
    Corners[2][0] = 0.5;
    Corners[2][1] = 0.866025;

    /* define barycentric coordinates */
    Scalar Barycentric_Coord[2];
    Barycentric_Coord[0] = Corners[0][0] * c1c + Corners[1][0] * c2c + Corners[2][0] * c3c;
    Barycentric_Coord[1] = Corners[0][1] * c1c + Corners[1][1] * c2c + Corners[2][1] * c3c;

    /* component 1C, 2C, 3C, converted to index of the "corners" */
    Scalar New_Coord[2];
    New_Coord[0] = Corners[uq_eigval_comp-1][0];
    New_Coord[1] = Corners[uq_eigval_comp-1][1];

    /* calculate perturbed barycentric coordinates */
    Barycentric_Coord[0] = Barycentric_Coord[0] + (uq_delta_b) * (New_Coord[0] - Barycentric_Coord[0]);
    Barycentric_Coord[1] = Barycentric_Coord[1] + (uq_delta_b) * (New_Coord[1] - Barycentric_Coord[1]);

    /* rebuild c1c,c2c,c3c based on perturbed barycentric coordinates */
    c3c = Barycentric_Coord[1] / Corners[2][1];
    c1c = Barycentric_Coord[0] - Corners[2][0] * c3c;
    c2c = 1 - c1c - c3c;

    /* build new anisotropy eigenvalues */
    Eig_Val[0] = (c3c - 1) / 3.0;
    Eig_Val[1] = 0.5 *c2c + Eig_Val[0];
    Eig_Val[2] = c1c + Eig_Val[1];

    /* permute eigenvectors if required */
    if (uq_permute) {
      for (size_t jDim = 0; jDim < 3; jDim++)
        swap(Eig_Vec[0][jDim], Eig_Vec[2][jDim]);
    }

    CBlasStructure::EigenRecomposition(A_ij, Eig_Vec, Eig_Val, 3);

    /* compute perturbed Reynolds stress matrix; using under-relaxation factor (uq_urlx)*/
    for (size_t iDim = 0; iDim < 3; iDim++) {
      for (size_t jDim = 0; jDim < 3; jDim++) {
        auto delta_ij = (jDim==iDim)? 1.0 : 0.0;
        MeanPerturbedRSM[iDim][jDim] = 2.0 * turb_ke * (A_ij[iDim][jDim] + 1.0/3.0 * delta_ij);
        MeanPerturbedRSM[iDim][jDim] = MeanReynoldsStress[iDim][jDim] +
          uq_urlx*(MeanPerturbedRSM[iDim][jDim] - MeanReynoldsStress[iDim][jDim]);
      }
    }
  }

  /*!
   * \brief Project average gradient onto normal (with or w/o correction) for viscous fluxes of scalar quantities.
   * \param[in] nDim - Dimension of the flow problem, 2 or 3.
   * \param[in] nVar - Number of variables.
   * \param[in] normal - Area vector.
   * \param[in] coord_i - Coordinate of point i.
   * \param[in] coord_j - Coordinate of point j.
   * \param[in] grad_i - Gradient at point i.
   * \param[in] grad_j - Gradient at point j.
   * \param[in] correct - Correct
   * \param[in] var_i - Variable at point i.
   * \param[in] var_j - Variable at point j.
   * \param[out] projNormal - Average gradient projected on normal.
   * \param[out] projCorrected - Corrected gradient projected on normal.
   * \return (Edge_Vector DOT normal) / |Edge_Vector|^2.
   */
  template<class Vec1, class Vec2, class Mat>
  FORCEINLINE static su2double ComputeProjectedGradient(int nDim, int nVar, const Vec1& normal,
                                                        const Vec1& coord_i, const Vec1& coord_j,
                                                        const Mat& grad_i, const Mat& grad_j,
                                                        bool correct,
                                                        const Vec2& var_i, const Vec2& var_j,
                                                        su2double* projNormal,
                                                        su2double* projCorrected) {
    assert(nDim == 2 || nDim == 3);
    nDim = (nDim > 2)? 3 : 2;
    su2double edgeVec[MAXNDIM], dist_ij_2 = 0.0, proj_vector_ij = 0.0;

    for (int iDim = 0; iDim < nDim; iDim++) {
      edgeVec[iDim] = coord_j[iDim] - coord_i[iDim];
      dist_ij_2 += pow(edgeVec[iDim], 2);
      proj_vector_ij += edgeVec[iDim] * normal[iDim];
    }
    proj_vector_ij /= max(dist_ij_2,EPS);

    /*--- Mean gradient approximation. ---*/
    for (int iVar = 0; iVar < nVar; iVar++) {
      projNormal[iVar] = 0.0;
      su2double edgeProj = 0.0;

      for (int iDim = 0; iDim < nDim; iDim++) {
        su2double meanGrad = 0.5 * (grad_i[iVar][iDim] + grad_j[iVar][iDim]);
        projNormal[iVar] += meanGrad * normal[iDim];
        if (correct) edgeProj += meanGrad * edgeVec[iDim];
      }

      projCorrected[iVar] = projNormal[iVar];
      if (correct) projCorrected[iVar] -= (edgeProj - (var_j[iVar]-var_i[iVar])) * proj_vector_ij;
    }

    return proj_vector_ij;
  }

  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first Menter blending function at point i.
   * \param[in] val_F1_j - Value of the first Menter blending function at point j.
   */
  virtual void SetF1blending(su2double val_F1_i, su2double val_F1_j) {/* empty */};

  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F2_i - Value of the second Menter blending function at point i.
   */
  virtual void SetF2blending(su2double val_F2_i) {/* empty */};

  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   */
  virtual void SetCrossDiff(su2double val_CDkw_i) {/* empty */};

  /*!
   * \brief Set the value of the effective intermittency for the LM model.
   * \param[in] intermittency_eff_i - Value of the effective intermittency at point i.
   */
  void SetIntermittencyEff(su2double val_intermittency_eff_i) {
    intermittency_eff_i = val_intermittency_eff_i;
  }

  /*!
   * \brief Set the value of the intermittency for the LM model.
   * \param[in] intermittency_i - Value of the intermittency at point i.
   */
  void SetIntermittency(su2double val_intermittency_i) {
    intermittency_i = val_intermittency_i;
  }

  /*!
   * \brief Get the value of the effective intermittency for the transition model.
   * \param[in] intermittency_eff_i - Value of the effective intermittency at point i.
   */
  su2double GetIntermittencyEff() const { return intermittency_eff_i; }

  /*!
   * \brief Set the gradient of the auxiliary variables.
   * \param[in] val_auxvar_grad_i - Gradient of the auxiliary variable at point i.
   * \param[in] val_auxvar_grad_j - Gradient of the auxiliary variable at point j.
   */
  inline void SetAuxVarGrad(CMatrixView<const su2double> val_auxvar_grad_i,
                            CMatrixView<const su2double> val_auxvar_grad_j) {
    AuxVar_Grad_i = val_auxvar_grad_i;
    AuxVar_Grad_j = val_auxvar_grad_j;
  }

  /*!
   * \brief Set the diffusion coefficient
   * \param[in] val_diffusioncoeff_i - Value of the diffusion coefficients at i.
   * \param[in] val_diffusioncoeff_j - Value of the diffusion coefficients at j
   */
  inline void SetDiffusionCoeff(const su2double* val_diffusioncoeff_i,
                                const su2double* val_diffusioncoeff_j) {
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
   * \brief Set the thermal conductivity (translational/rotational)
   * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
   * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  inline void SetThermalConductivity_ve(su2double val_thermal_conductivity_ve_i,
                                     su2double val_thermal_conductivity_ve_j) {
    Thermal_Conductivity_ve_i = val_thermal_conductivity_ve_i;
    Thermal_Conductivity_ve_j = val_thermal_conductivity_ve_j;
  }

  /*!
   * \brief Set the specifc heat c_p.
   * \param[in] val_specific_heat_i - Value of the specific heat at point i.
   * \param[in] val_specific_heat_j - Value of the specific heat at point j.
   */
  inline void SetSpecificHeat(su2double val_specific_heat_i,
                              su2double val_specific_heat_j) {
    Cp_i = val_specific_heat_i;
    Cp_j = val_specific_heat_j;
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
   * \brief Set the value of the roughness from the nearest wall.
   * \param[in] val_dist_i - Value of of the roughness of the nearest wall from point i
   * \param[in] val_dist_j - Value of of the roughness of the nearest wall from point j
   */
  void SetRoughness(su2double val_roughness_i, su2double val_roughness_j) {
    roughness_i = val_roughness_i;
    roughness_j = val_roughness_j;
  }

  /*!
   * \brief Set coordinates of the points.
   * \param[in] val_coord_i - Coordinates of the point i.
   * \param[in] val_coord_j - Coordinates of the point j.
   */
  inline void SetCoord(const su2double *val_coord_i, const su2double *val_coord_j) {
    Coord_i = val_coord_i;
    Coord_j = val_coord_j;
  }

  /*!
   * \brief Set the velocity of the computational grid.
   * \param[in] val_gridvel_i - Grid velocity of the point i.
   * \param[in] val_gridvel_j - Grid velocity of the point j.
   */
  inline void SetGridVel(const su2double *val_gridvel_i, const su2double *val_gridvel_j) {
    GridVel_i = val_gridvel_i;
    GridVel_j = val_gridvel_j;
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
  inline void SetUndivided_Laplacian(const su2double *val_und_lapl_i, const su2double *val_und_lapl_j) {
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
  inline void SetNormal(const su2double *val_normal) { Normal = val_normal; }

  /*!
   * \brief Set the value of the volume of the control volume.
   * \param[in] val_volume Volume of the control volume.
   */
  inline void SetVolume(su2double val_volume) { Volume = val_volume; }

  /*!
   * \brief Set the value of AvgVolume Variable.
   * \param[in] val_avgvolume AvgVolume Variable.
   */
  inline void SetAvgVolume(su2double val_avgvolume) { AvgVolume = val_avgvolume; }

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
   * \brief Compute the projected inviscid flux vector.
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_enthalpy - Pointer to the enthalpy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidProjFlux(const su2double *val_density, const su2double *val_velocity,
                           const su2double *val_pressure, const su2double *val_enthalpy,
                           const su2double *val_normal, su2double *val_Proj_Flux) const;

  /*!
   * \brief Compute the projected inviscid flux vector for incompresible simulations
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidIncProjFlux(const su2double *val_density, const su2double *val_velocity,
                              const su2double *val_pressure, const su2double *val_betainc2,
                              const su2double *val_enthalpy, const su2double *val_normal,
                              su2double *val_Proj_Flux) const;

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices.
   * \param[in] val_velocity Pointer to the velocity.
   * \param[in] val_energy Value of the energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(const su2double *val_velocity, const su2double *val_energy,
                          const su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_tensor) const;

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices (incompressible).
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidIncProjJac(const su2double *val_density, const su2double *val_velocity,
                             const su2double *val_betainc2, const su2double *val_normal,
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
  void GetInviscidIncProjJac(const su2double *val_density,
                             const su2double *val_velocity,
                             const su2double *val_betainc2,
                             const su2double *val_cp,
                             const su2double *val_temperature,
                             const su2double *val_dRhodT,
                             const su2double *val_normal,
                             su2double val_scale,
                             su2double **val_Proj_Jac_Tensor) const;

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
  void GetPreconditioner(const su2double *val_density,
                         const su2double *val_velocity,
                         const su2double *val_betainc2,
                         const su2double *val_cp,
                         const su2double *val_temperature,
                         const su2double *val_drhodt,
                         su2double **val_Precon) const;

  /*!
   * \brief Compute the projection of the preconditioned inviscid Jacobian matrices.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetPreconditionedProjJac(const su2double *val_density,
                                const su2double *val_velocity,
                                const su2double *val_betainc2,
                                const su2double *val_normal,
                                su2double **val_Proj_Jac_Tensor) const;

  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices for general fluid model.
   * \param[in] val_velocity Pointer to the velocity.
   * \param[in] val_energy Value of the energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(const su2double *val_velocity, const su2double *val_enthalphy,
                          const su2double *val_chi, const su2double *val_kappa,
                          const su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_tensor) const;

  /*!
   * \brief Mapping between primitives variables P and conservatives variables C.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_PrimVar - Mean Value of the secondary variables.
   * \param[out] val_Jac_PC - Pointer to the Jacobian dPdC.
   */
  void GetPrimitive2Conservative (const su2double *val_Mean_PrimVar,
                                  const su2double *val_Mean_SecVar,
                                  su2double **val_Jac_PC) const;

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
  void GetPMatrix(const su2double *val_density, const su2double *val_velocity,
                  const su2double *val_soundspeed, const su2double *val_enthalpy,
                  const su2double *val_chi, const su2double *val_kappa,
                  const su2double *val_normal, su2double **val_p_tensor) const;

  /*!
   * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in
   *        the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPMatrix(const su2double *val_density, const su2double *val_velocity,
                  const su2double *val_soundspeed, const su2double *val_normal,
                  su2double **val_p_tensor) const;

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
                    su2double val_density, const su2double* val_velocity,
                    su2double** val_invR_invPe) const;

  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_pressure - value of the pressure.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[in] val_velocity - value of the velocity.
   * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetRMatrix(su2double val_pressure, su2double val_soundspeed,
                  su2double val_density, const su2double* val_velocity,
                  su2double** val_invR_invPe) const;
  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[out] R_Matrix - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetRMatrix(su2double val_soundspeed, su2double val_density, su2double **R_Matrix) const;

  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[out] L_Matrix - Pointer to the matrix of conversion from conserved to entropic variables.
   */
  void GetLMatrix(su2double val_soundspeed, su2double val_density, su2double **L_Matrix) const;

  /*!
   * \brief Computation of the flow Residual Jacobian Matrix for Non Reflecting BC.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[out] R_c - Residual Jacoboan Matrix
   * \param[out] R_c_inv- inverse of the Residual Jacobian Matrix .
   */
  void ComputeResJacobianGiles(CFluidModel *FluidModel, su2double pressure, su2double density, const su2double *turboVel,
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
  void GetCharJump(su2double val_soundspeed, su2double val_density, const su2double *prim_jump, su2double *char_jump) const;

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
                          su2double rB2a2, const su2double* val_Lambda, const su2double* val_normal, su2double** val_absPeJac) const;

  /*!
   * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
   * in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(su2double **val_invp_tensor, const su2double *val_density,
                      const su2double *val_velocity, const su2double *val_soundspeed,
                      const su2double *val_chi, const su2double *val_kappa,
                      const su2double *val_normal) const;

  /*!
   * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(const su2double *val_density, const su2double *val_velocity,
                      const su2double *val_soundspeed, const su2double *val_normal,
                      su2double **val_invp_tensor) const;

  /*!
   * \brief Compute viscous residual and jacobian.
   */
  void GetAdjViscousFlux_Jac(su2double Pressure_i, su2double Pressure_j, su2double Density_i, su2double Density_j,
                             su2double ViscDens_i, su2double ViscDens_j, const su2double *Velocity_i, const su2double *Velocity_j,
                             su2double sq_vel_i, su2double sq_vel_j,
                             su2double XiDens_i, su2double XiDens_j, su2double **Mean_GradPhi, const su2double *Mean_GradPsiE,
                             su2double dPhiE_dn, const su2double *Normal, const su2double *Edge_Vector, su2double dist_ij_2,
                             su2double *val_residual_i, su2double *val_residual_j,
                             su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                             su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, bool implicit) const;

  /*!
   * \brief Computation of the projected inviscid lambda (eingenvalues).
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_Lambda_Vector - Pointer to Lambda matrix.
   */
  void GetJacInviscidLambda_fabs(const su2double *val_velocity, su2double val_soundspeed,
                                 const su2double *val_normal, su2double *val_Lambda_Vector) const;

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
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ComputeResidual(su2double *val_resconv, su2double *val_resvisc,
                               su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                               CConfig *config) {}

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
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual ResidualType<> ComputeAxisymmetric(const CConfig* config) { return ResidualType<>(nullptr,nullptr,nullptr); }

  /*!
   * \overload For numerics classes that store the residual/flux and Jacobians internally.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  inline virtual ResidualType<> ComputeVibRelaxation(const CConfig* config) { return ResidualType<>(nullptr,nullptr,nullptr); }

  /*!
   * \brief Calculation of the chemistry source term
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  inline virtual ResidualType<> ComputeChemistry(const CConfig* config) { return ResidualType<>(nullptr,nullptr,nullptr); }

  /*!
   * \brief Check if residual constains a NaN value
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_residual - residual of the numeric function.
   * \param[out] ERR - Presencse of NaN in vector
   */
  static bool CheckResidualNaNs(bool implicit, int nVar, const ResidualType<> residual) {

    bool ERR = false;
    const bool jac_j = residual.jacobian_j != nullptr;

    for (auto iVar = 0; iVar<nVar; iVar++){
      if (std::isnan(SU2_TYPE::GetValue(residual[iVar]))) ERR = true;

      if (implicit) {
        for (auto jVar = 0; jVar < nVar; jVar++){
          if (std::isnan(SU2_TYPE::GetValue(residual.jacobian_i[iVar][jVar]))) ERR = true;
          if ((jac_j) && (std::isnan(SU2_TYPE::GetValue(residual.jacobian_j[iVar][jVar])))) ERR = true;
        }
      }
    }
    return ERR;
  }

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
  inline virtual su2double Compute_Averaged_NodalStress(CElement *element_container, const CConfig* config) { return 0; }

  /*!
   * \brief Computes a basis of orthogonal vectors from a supplied vector
   * \param[in] config - Normal vector
   */
  void CreateBasis(const su2double *val_Normal, su2double* l, su2double* m);

  /*!
   * \brief Set the value of the Tauwall
   * \param[in] val_tauwall_i - Tauwall at point i
   * \param[in] val_tauwall_j - Tauwall at point j
   */
  inline virtual void SetTau_Wall(su2double val_tauwall_i, su2double val_tauwall_j) { }

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
  inline virtual void SetRadVarGradient(CMatrixView<const su2double> val_radvar_grad_i,
                                        CMatrixView<const su2double> val_radvar_grad_j) { }

  /*!
   * \brief Set the gradient of the radiation variables.
   * \param[in] val_radvar_source - Source term (and jacobian term) of the radiative heat transfer.
   */
  inline void SetRadVarSource(const su2double *val_radvar_source) { RadVar_Source = val_radvar_source; }

  /*!
   * \brief Set the pressure derivatives.
   * \param[in] val_dPdU_i - pressure derivatives at i.
   * \param[in] val_dPdU_j - pressure derivatives at j.
   */
  virtual inline void SetdPdU(su2double *val_dPdU_i, su2double *val_dPdU_j)       { }

  /*!
   * \brief Set the temperature derivatives.
   * \param[in] val_dTdU_i - temperature derivatives at i.
   * \param[in] val_dTdU_j - temperature derivatives at j.
   */
  virtual inline void SetdTdU(su2double *val_dTdU_i, su2double *val_dTdU_j)       { }

  /*!
   * \brief Set the vib-el temperture derivatives.
   * \param[in] val_dTvedU_i - t_ve derivatives at i.
   * \param[in] val_dTvedU_j - t_ve derivatives at j.
   */
  virtual inline void SetdTvedU(su2double *val_dTvedU_i, su2double *val_dTvedU_j) { }

  /*!
   * \brief Set the vib-elec energy.
   * \param[in] val_Eve_i - vib-el energy at i.
   * \param[in] val_Eve_j - vib-el energy at j.
   */
  virtual inline void SetEve(su2double *val_Eve_i, su2double *val_Eve_j)          { }

  /*!
   * \brief Set the vib-elec specific heat.
   * \param[in] val_Cvve_i - Cvve at i.
   * \param[in] val_Cvve_j - Cvve at j.
   */
  virtual inline void SetCvve(su2double *val_Cvve_i, su2double *val_Cvve_j)       { }

  /*!
   * \brief Set the ratio of specific heats.
   * \param[in] val_Gamma_i - Gamma at i.
   * \param[in] val_Gamma_j - Gamma at j.
   */
  virtual inline void SetGamma(su2double val_Gamma_i, su2double val_Gamma_j)       { }

  /*!
   * \brief Set massflow, heatflow & inlet temperature for streamwise periodic flow.
   * \param[in] SolverSPvals - Struct holding the values.
   */
  virtual void SetStreamwisePeriodicValues(const StreamwisePeriodicValues SolverSPvals) { }

  /*!
   * \brief SetMassFlux
   * \param[in] val_MassFlux: Mass flux across the edge
   */
  inline void SetMassFlux(const su2double val_MassFlux) { MassFlux = val_MassFlux; }

  /*!
   * \brief Obtain information on bounded scalar problem
   * \return is_bounded_scalar : scalar solver uses bounded scalar convective transport
   */
  inline bool GetBoundedScalar() const { return bounded_scalar;}
};

/*!
 * /brief Alias for a "do nothing" source term class
 */
using CSourceNothing = CNumerics;
