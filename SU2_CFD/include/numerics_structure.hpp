/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/gauss_structure.hpp"
#include "../../Common/include/element_structure.hpp"
#include "variable_structure.hpp"

using namespace std;

/*!
 * \class CNumerics
 * \brief Class for defining the numerical methods.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 */
class CNumerics {
protected:
  unsigned short nDim, nVar;  /*!< \brief Number of dimensions and variables. */
  su2double Gamma;        /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;    /*!< \brief Fluids's Gamma - 1.0  . */
  su2double Gas_Constant;         /*!< \brief Gas constant. */
  su2double *Vector; /*!< \brief Auxiliary vector. */
  su2double *Enthalpy_formation;
  su2double Prandtl_Lam;        /*!< \brief Laminar Prandtl's number. */
  su2double Prandtl_Turb;    /*!< \brief Turbulent Prandtl's number. */
  
public:
  
  su2double
  **Flux_Tensor,  /*!< \brief Flux tensor (used for viscous and inviscid purposes. */
  *Proj_Flux_Tensor;    /*!< \brief Flux tensor projected in a direction. */
  
  su2double
  **tau,    /*!< \brief Viscous stress tensor. */
  **delta;      /*!< \brief Identity matrix. */
  su2double
  *Diffusion_Coeff_i, /*!< \brief Species diffusion coefficients at point i. */
  *Diffusion_Coeff_j; /*!< \brief Species diffusion coefficients at point j. */
  su2double Laminar_Viscosity_i,  /*!< \brief Laminar viscosity at point i. */
  Laminar_Viscosity_j,    /*!< \brief Laminar viscosity at point j. */
  Laminar_Viscosity_id,  /*!< \brief Variation of laminar viscosity at point i. */
  Laminar_Viscosity_jd;    /*!< \brief Variation of laminar viscosity at point j. */
  su2double Thermal_Conductivity_i, /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_j, /*!< \brief Thermal conductivity at point j. */
  Thermal_Conductivity_ve_i, /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_ve_j; /*!< \brief Thermal conductivity at point j. */
  su2double Cp_i, /*!< \brief Cp at point i. */
  Cp_j;         /*!< \brief Cp at point j. */
  su2double *Theta_v; /*!< \brief Characteristic vibrational temperature */
  su2double Eddy_Viscosity_i,  /*!< \brief Eddy viscosity at point i. */
  Eddy_Viscosity_j;      /*!< \brief Eddy viscosity at point j. */
  su2double turb_ke_i,  /*!< \brief Turbulent kinetic energy at point i. */
  turb_ke_j;      /*!< \brief Turbulent kinetic energy at point j. */
  su2double Pressure_i,  /*!< \brief Pressure at point i. */
  Pressure_j;      /*!< \brief Pressure at point j. */
  su2double GravityForce_i,  /*!< \brief Gravity force at point i. */
  GravityForce_j;      /*!< \brief Gravity force at point j. */
  su2double Density_i,  /*!< \brief Density at point i. */
  Density_j;      /*!< \brief Density at point j. */
  su2double DensityInc_i,  /*!< \brief Incompressible density at point i. */
  DensityInc_j;      /*!< \brief Incompressible density at point j. */
  su2double BetaInc2_i,  /*!< \brief Beta incompressible at point i. */
  BetaInc2_j;      /*!< \brief Beta incompressible at point j. */
  su2double Lambda_i,  /*!< \brief Spectral radius at point i. */
  Lambda_j;      /*!< \brief Spectral radius at point j. */
  su2double LambdaComb_i,  /*!< \brief Spectral radius at point i. */
  LambdaComb_j;      /*!< \brief Spectral radius at point j. */
  su2double SoundSpeed_i,  /*!< \brief Sound speed at point i. */
  SoundSpeed_j;      /*!< \brief Sound speed at point j. */
  su2double Enthalpy_i,  /*!< \brief Enthalpy at point i. */
  Enthalpy_j;      /*!< \brief Enthalpy at point j. */
  su2double dist_i,  /*!< \brief Distance of point i to the nearest wall. */
  dist_j;      /*!< \brief Distance of point j to the nearest wall. */
  su2double Temp_i,  /*!< \brief Temperature at point i. */
  Temp_j;      /*!< \brief Temperature at point j. */
  su2double *Temp_tr_i, /*!< \brief Temperature transl-rot at point i. */
  *Temp_tr_j;/*!< \brief Temperature transl-rot at point j. */
  su2double *Temp_vib_i, /*!< \brief Temperature vibrational at point i. */
  *Temp_vib_j;/*!< \brief Temperature vibrational at point j. */
  su2double *Und_Lapl_i, /*!< \brief Undivided laplacians at point i. */
  *Und_Lapl_j;    /*!< \brief Undivided laplacians at point j. */
  su2double Sensor_i,  /*!< \brief Pressure sensor at point i. */
  Sensor_j;      /*!< \brief Pressure sensor at point j. */
  su2double *GridVel_i,  /*!< \brief Grid velocity at point i. */
  *GridVel_j;      /*!< \brief Grid velocity at point j. */
  su2double *U_i,    /*!< \brief Vector of conservative variables at point i. */
  *U_id,    /*!< \brief Vector of derivative of conservative variables at point i. */
  *UZeroOrder_i,  /*!< \brief Vector of conservative variables at point i without reconstruction. */
  *U_j,        /*!< \brief Vector of conservative variables at point j. */
  *UZeroOrder_j,  /*!< \brief Vector of conservative variables at point j without reconstruction. */
  *U_jd,        /*!< \brief Vector of derivative of conservative variables at point j. */
  *U_0,        /*!< \brief Vector of conservative variables at node 0. */
  *U_1,        /*!< \brief Vector of conservative variables at node 1. */
  *U_2,        /*!< \brief Vector of conservative variables at node 2. */
  *U_3;        /*!< \brief Vector of conservative variables at node 3. */
  su2double *V_i,    /*!< \brief Vector of primitive variables at point i. */
  *V_j;        /*!< \brief Vector of primitive variables at point j. */
  su2double *S_i,    /*!< \brief Vector of secondary variables at point i. */
  *S_j;        /*!< \brief Vector of secondary variables at point j. */
  su2double *Psi_i,    /*!< \brief Vector of adjoint variables at point i. */
  *Psi_j;        /*!< \brief Vector of adjoint variables at point j. */
  su2double *DeltaU_i,  /*!< \brief Vector of linearized variables at point i. */
  *DeltaU_j;      /*!< \brief Vector of linearized variables at point j. */
  su2double *TurbVar_i,  /*!< \brief Vector of turbulent variables at point i. */
  *TurbVar_id,  /*!< \brief Vector of derivative of turbulent variables at point i. */
  *TurbVar_j,      /*!< \brief Vector of turbulent variables at point j. */
  *TurbVar_jd;  /*!< \brief Vector of derivative of turbulent variables at point j. */
  su2double *TransVar_i,  /*!< \brief Vector of turbulent variables at point i. */
  *TransVar_j;      /*!< \brief Vector of turbulent variables at point j. */
  su2double *TurbPsi_i,  /*!< \brief Vector of adjoint turbulent variables at point i. */
  *TurbPsi_j;      /*!< \brief Vector of adjoint turbulent variables at point j. */
  su2double **ConsVar_Grad_i,  /*!< \brief Gradient of conservative variables at point i. */
  **ConsVar_Grad_j,      /*!< \brief Gradient of conservative variables at point j. */
  **ConsVar_Grad_0,      /*!< \brief Gradient of conservative variables at point 0. */
  **ConsVar_Grad_1,      /*!< \brief Gradient of conservative variables at point 1. */
  **ConsVar_Grad_2,      /*!< \brief Gradient of conservative variables at point 2. */
  **ConsVar_Grad_3,      /*!< \brief Gradient of conservative variables at point 3. */
  **ConsVar_Grad;        /*!< \brief Gradient of conservative variables which is a scalar. */
  su2double **PrimVar_Grad_i,  /*!< \brief Gradient of primitive variables at point i. */
  **PrimVar_Grad_j;      /*!< \brief Gradient of primitive variables at point j. */
  su2double **PsiVar_Grad_i,    /*!< \brief Gradient of adjoint variables at point i. */
  **PsiVar_Grad_j;      /*!< \brief Gradient of adjoint variables at point j. */
  su2double **TurbVar_Grad_i,  /*!< \brief Gradient of turbulent variables at point i. */
  **TurbVar_Grad_j;      /*!< \brief Gradient of turbulent variables at point j. */
  su2double **TransVar_Grad_i,  /*!< \brief Gradient of turbulent variables at point i. */
  **TransVar_Grad_j;      /*!< \brief Gradient of turbulent variables at point j. */
  su2double **TurbPsi_Grad_i,  /*!< \brief Gradient of adjoint turbulent variables at point i. */
  **TurbPsi_Grad_j;      /*!< \brief Gradient of adjoint turbulent variables at point j. */
  su2double *AuxVar_Grad_i,    /*!< \brief Gradient of an auxiliary variable at point i. */
  *AuxVar_Grad_j;        /*!< \brief Gradient of an auxiliary variable at point i. */
  su2double *Coord_i,  /*!< \brief Cartesians coordinates of point i. */
  *Coord_j,      /*!< \brief Cartesians coordinates of point j. */
  *Coord_0,      /*!< \brief Cartesians coordinates of point 0 (Galerkin method, triangle). */
  *Coord_1,      /*!< \brief Cartesians coordinates of point 1 (Galerkin method, tetrahedra). */
  *Coord_2,      /*!< \brief Cartesians coordinates of point 2 (Galerkin method, triangle). */
  *Coord_3;      /*!< \brief Cartesians coordinates of point 3 (Galerkin method, tetrahedra). */
  unsigned short Neighbor_i,  /*!< \brief Number of neighbors of the point i. */
  Neighbor_j;          /*!< \brief Number of neighbors of the point j. */
  su2double *Normal,  /*!< \brief Normal vector, it norm is the area of the face. */
  *UnitNormal,    /*!< \brief Unitary normal vector. */
  *UnitNormald;    /*!< \brief derivatve of unitary normal vector. */
  su2double TimeStep,    /*!< \brief Time step useful in dual time method. */
  Area,        /*!< \brief Area of the face i-j. */
  Volume;        /*!< \brief Volume of the control volume around point i. */
  su2double Volume_n,  /*!< \brief Volume of the control volume at time n. */
  Volume_nM1,    /*!< \brief Volume of the control volume at time n-1. */
  Volume_nP1;    /*!< \brief Volume of the control volume at time n+1. */
  su2double *U_n,  /*!< \brief Vector of conservative variables at time n. */
  *U_nM1,    /*!< \brief Vector of conservative variables at time n-1. */
  *U_nP1;    /*!< \brief Vector of conservative variables at time n+1. */
  su2double vel2_inf; /*!< \brief value of the square of freestream speed. */
  su2double *WindGust_i,  /*!< \brief Wind gust at point i. */
  *WindGust_j;      /*!< \brief Wind gust at point j. */
  su2double *WindGustDer_i,  /*!< \brief Wind gust derivatives at point i. */
  *WindGustDer_j;      /*!< \brief Wind gust derivatives at point j. */
  su2double *Vorticity_i, *Vorticity_j;  /*!< \brief Vorticity. */
  su2double StrainMag_i, StrainMag_j;   /*!< \brief Strain rate magnitude. */
  
  su2double *l, *m;

  su2double **MeanReynoldsStress;
  su2double **MeanPerturbedRSM;
  bool using_uq;
  su2double PerturbedStrainMag;
  
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
  CNumerics(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CNumerics(void);
  
  /*!
   * \brief Compute the determinant of a 3 by 3 matrix.
   * \param[in] val_matrix 3 by 3 matrix.
   * \result Determinant of the matrix
   */
  su2double Determinant_3x3(su2double A00, su2double A01, su2double A02,
                            su2double A10, su2double A11, su2double A12,
                            su2double A20, su2double A21, su2double A22);
  
  /*!
   * \brief Set the solution at different times.
   * \param[in] val_u_nM1 Conservative solution at time n-1.
   * \param[in] val_u_n Conservative solution at time n.
   * \param[in] val_u_nP1 Conservative solution at time n+1.
   */
  void SetPastSol(su2double *val_u_nM1, su2double *val_u_n, su2double *val_u_nP1);
  
  /*!
   * \brief Set the control volume at different times.
   * \param[in] val_volume_nM1 - Control volume at time n-1.
   * \param[in] val_volume_n - Control volume at time n.
   * \param[in] val_volume_nP1 - Control volume at time n+1.
   */
  void SetPastVolume(su2double val_volume_nM1, su2double val_volume_n, su2double val_volume_nP1);
  
  /*!
   * \brief Set the time step.
   * \param[in] val_timestep - Value of the time step.
   */
  void SetTimeStep(su2double val_timestep);
  
  /*!
   * \brief Get the Preconditioning Beta.
   * \return val_Beta - Value of the low Mach Preconditioner.
   */
  virtual su2double GetPrecond_Beta();
  
  /*!
   * \brief Set the freestream velocity square.
   * \param[in] SetVelocity2_Inf - Value of the square of the freestream velocity.
   */
  void SetVelocity2_Inf(su2double val_velocity2);
  
  /*!
   * \brief Set the value of the vorticity
   * \param[in] val_vorticity - Value of the vorticity.
   */
  void SetVorticity(su2double *val_vorticity_i, su2double *val_vorticity_j);
  
  /*!
   * \brief Set the value of the rate of strain magnitude.
   * \param[in] val_StrainMag_i - Value of the magnitude of rate of strain at point i.
   * \param[in] val_StrainMag_j - Value of the magnitude of rate of strain at point j.
   */
  void SetStrainMag(su2double val_strainmag_i, su2double val_strainmag_j);
  
  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_i - Value of the conservative variable at point i.
   * \param[in] val_u_j - Value of the conservative variable at point j.
   */
  void SetConservative(su2double *val_u_i, su2double *val_u_j);
  
  /*!
   * \brief Set the value of the conservative variables withour reconstruction.
   * \param[in] val_u_i - Value of the conservative variable at point i.
   * \param[in] val_u_j - Value of the conservative variable at point j.
   */
  void SetConservative_ZeroOrder(su2double *val_u_i, su2double *val_u_j);
  
  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  void SetPrimitive(su2double *val_v_i, su2double *val_v_j);
  
  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  void SetSecondary(su2double *val_s_i, su2double *val_s_j);
  
  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_0 - Value of the conservative variable at point 0.
   * \param[in] val_u_1 - Value of the conservative variable at point 1.
   * \param[in] val_u_2 - Value of the conservative variable at point 2.
   */
  void SetConservative(su2double *val_u_0, su2double *val_u_1, su2double *val_u_2);
  
  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_0 - Value of the conservative variable at point 0.
   * \param[in] val_u_1 - Value of the conservative variable at point 1.
   * \param[in] val_u_2 - Value of the conservative variable at point 2.
   * \param[in] val_u_3 - Value of the conservative variable at point 3.
   */
  void SetConservative(su2double *val_u_0, su2double *val_u_1, su2double *val_u_2, su2double *val_u_3);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_i - Gradient of the conservative variable at point i.
   * \param[in] val_consvar_grad_j - Gradient of the conservative variable at point j.
   */
  void SetConsVarGradient(su2double **val_consvar_grad_i, su2double **val_consvar_grad_j);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
   * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
   * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
   */
  void SetConsVarGradient(su2double **val_consvar_grad_0,
                          su2double **val_consvar_grad_1,
                          su2double **val_consvar_grad_2);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
   * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
   * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
   * \param[in] val_consvar_grad_3 - Gradient of the conservative variable at point 3.
   */
  void SetConsVarGradient(su2double **val_consvar_grad_0,
                          su2double **val_consvar_grad_1,
                          su2double **val_consvar_grad_2,
                          su2double **val_consvar_grad_3);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad - Gradient of the conservative variable which is a scalar.
   */
  void SetConsVarGradient(su2double **val_consvar_grad);
  
  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
   * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
   */
  void SetPrimVarGradient(su2double **val_primvar_grad_i,
                          su2double **val_primvar_grad_j);
  
  /*!
   * \brief Set the value of the adjoint variable.
   * \param[in] val_psi_i - Value of the adjoint variable at point i.
   * \param[in] val_psi_j - Value of the adjoint variable at point j.
   */
  void SetAdjointVar(su2double *val_psi_i, su2double *val_psi_j);
  
  /*!
   * \brief Set the gradient of the adjoint variables.
   * \param[in] val_psivar_grad_i - Gradient of the adjoint variable at point i.
   * \param[in] val_psivar_grad_j - Gradient of the adjoint variable at point j.
   */
  void SetAdjointVarGradient(su2double **val_psivar_grad_i, su2double **val_psivar_grad_j);
  
  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
   * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
   */
  void SetTurbVar(su2double *val_turbvar_i, su2double *val_turbvar_j);
  
  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_transvar_i - Value of the turbulent variable at point i.
   * \param[in] val_transvar_j - Value of the turbulent variable at point j.
   */
  void SetTransVar(su2double *val_transvar_i, su2double *val_transvar_j);
  
  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  void SetTurbVarGradient(su2double **val_turbvar_grad_i, su2double **val_turbvar_grad_j);
  
  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  void SetTransVarGradient(su2double **val_transvar_grad_i, su2double **val_transvar_grad_j);

  /*!
   * \brief Set the value of the adjoint turbulent variable.
   * \param[in] val_turbpsivar_i - Value of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_j - Value of the adjoint turbulent variable at point j.
   */
  void SetTurbAdjointVar(su2double *val_turbpsivar_i, su2double *val_turbpsivar_j);
  
  /*!
   * \brief Set the gradient of the adjoint turbulent variables.
   * \param[in] val_turbpsivar_grad_i - Gradient of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_grad_j - Gradient of the adjoint turbulent variable at point j.
   */
  void SetTurbAdjointGradient (su2double **val_turbpsivar_grad_i, su2double **val_turbpsivar_grad_j);
  
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
  void SetAuxVarGrad(su2double *val_auxvargrad_i, su2double *val_auxvargrad_j);
  
  /*!
   * \brief Set the diffusion coefficient
   * \param[in] val_diffusioncoeff_i - Value of the diffusion coefficients at i.
   * \param[in] val_diffusioncoeff_j - Value of the diffusion coefficients at j
   */
  void SetDiffusionCoeff(su2double* val_diffusioncoeff_i,
                         su2double* val_diffusioncoeff_j);
  
  /*!
   * \brief Set the laminar viscosity.
   * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
   * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
   */
  void SetLaminarViscosity(su2double val_laminar_viscosity_i,
                           su2double val_laminar_viscosity_j);
  
  /*!
   * \brief Set the thermal conductivity (translational/rotational)
   * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
   * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  void SetThermalConductivity(su2double val_thermal_conductivity_i,
                              su2double val_thermal_conductivity_j);
  
  /*!
   * \brief Set the thermal conductivity (translational/rotational)
   * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
   * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  void SetThermalConductivity_ve(su2double val_thermal_conductivity_ve_i,
                                 su2double val_thermal_conductivity_ve_j);
  
  /*!
   * \brief Set the eddy viscosity.
   * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
   * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
   */
  void SetEddyViscosity(su2double val_eddy_viscosity_i,
                        su2double val_eddy_viscosity_j);
  
  /*!
   * \brief Set the turbulent kinetic energy.
   * \param[in] val_turb_ke_i - Value of the turbulent kinetic energy at point i.
   * \param[in] val_turb_ke_j - Value of the turbulent kinetic energy at point j.
   */
  void SetTurbKineticEnergy(su2double val_turb_ke_i, su2double val_turb_ke_j);
  
  /*!
   * \brief Set the value of the distance from the nearest wall.
   * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
   * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
   */
  void SetDistance(su2double val_dist_i, su2double val_dist_j);
  
  /*!
   * \brief Set coordinates of the points.
   * \param[in] val_coord_i - Coordinates of the point i.
   * \param[in] val_coord_j - Coordinates of the point j.
   */
  void SetCoord(su2double *val_coord_i, su2double *val_coord_j);
  
  /*!
   * \overload
   * \param[in] val_coord_0 - Coordinates of the point 0.
   * \param[in] val_coord_1 - Coordinates of the point 1.
   * \param[in] val_coord_2 - Coordinates of the point 2.
   */
  void SetCoord(su2double *val_coord_0, su2double *val_coord_1, su2double *val_coord_2);
  
  /*!
   * \overload
   * \param[in] val_coord_0 - Coordinates of the point 0.
   * \param[in] val_coord_1 - Coordinates of the point 1.
   * \param[in] val_coord_2 - Coordinates of the point 2.
   * \param[in] val_coord_3 - Coordinates of the point 3.
   */
  void SetCoord(su2double *val_coord_0, su2double *val_coord_1, su2double *val_coord_2,
                su2double *val_coord_3);
  
  /*!
   * \brief Set the velocity of the computational grid.
   * \param[in] val_gridvel_i - Grid velocity of the point i.
   * \param[in] val_gridvel_j - Grid velocity of the point j.
   */
  void SetGridVel(su2double *val_gridvel_i, su2double *val_gridvel_j);
  
  /*!
   * \brief Set the wind gust value.
   * \param[in] val_windgust_i - Wind gust of the point i.
   * \param[in] val_windgust_j - Wind gust of the point j.
   */
  void SetWindGust(su2double *val_windgust_i, su2double *val_windgust_j);
  
  /*!
   * \brief Set the wind gust derivatives values.
   * \param[in] val_windgust_i - Wind gust derivatives of the point i.
   * \param[in] val_windgust_j - Wind gust derivatives of the point j.
   */
  void SetWindGustDer(su2double *val_windgustder_i, su2double *val_windgustder_j);
  
  /*!
   * \brief Set the value of the pressure.
   * \param[in] val_pressure_i - Value of the pressure at point i.
   * \param[in] val_pressure_j - Value of the pressure at point j.
   */
  void SetPressure(su2double val_pressure_i, su2double val_pressure_j);
  
  /*!
   * \brief Set the value of the density for the incompressible solver.
   * \param[in] val_densityinc_i - Value of the pressure at point i.
   * \param[in] val_densityinc_j - Value of the pressure at point j.
   */
  void SetDensity(su2double val_densityinc_i, su2double val_densityinc_j);
  
  /*!
   * \brief Set the value of the beta for incompressible flows.
   * \param[in] val_betainc2_i - Value of beta for incompressible flows at point i.
   * \param[in] val_betainc2_j - Value of beta for incompressible flows at point j.
   */
  void SetBetaInc2(su2double val_betainc2_i, su2double val_betainc2_j);
  
  /*!
   * \brief Set the value of the sound speed.
   * \param[in] val_soundspeed_i - Value of the sound speed at point i.
   * \param[in] val_soundspeed_j - Value of the sound speed at point j.
   */
  void SetSoundSpeed(su2double val_soundspeed_i, su2double val_soundspeed_j);
  
  /*!
   * \brief Set the value of the temperature.
   * \param[in] val_temp_i - Value of the temperature at point i.
   * \param[in] val_temp_j - Value of the temperature at point j.
   */
  void SetTemperature(su2double val_temp_i, su2double val_temp_j);
  
  /*!
   * \brief Set the value of the enthalpy.
   * \param[in] val_enthalpy_i - Value of the enthalpy at point i.
   * \param[in] val_enthalpy_j - Value of the enthalpy at point j.
   */
  void SetEnthalpy(su2double val_enthalpy_i, su2double val_enthalpy_j);
  
  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] val_lambda_i - Value of the spectral radius at point i.
   * \param[in] val_lambda_j - Value of the spectral radius at point j.
   */
  void SetLambda(su2double val_lambda_i, su2double val_lambda_j);
  
  /*!
   * \brief Set the value of undivided laplacian.
   * \param[in] val_und_lapl_i Undivided laplacian at point i.
   * \param[in] val_und_lapl_j Undivided laplacian at point j.
   */
  void SetUndivided_Laplacian(su2double *val_und_lapl_i, su2double *val_und_lapl_j);
  
  /*!
   * \brief Set the value of the pressure sensor.
   * \param[in] val_sensor_i Pressure sensor at point i.
   * \param[in] val_sensor_j Pressure sensor at point j.
   */
  void SetSensor(su2double val_sensor_i, su2double val_sensor_j);
  
  /*!
   * \brief Set the number of neighbor to a point.
   * \param[in] val_neighbor_i - Number of neighbor to point i.
   * \param[in] val_neighbor_j - Number of neighbor to point j.
   */
  void SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j);
  
  /*!
   * \brief Set the value of the normal vector to the face between two points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void SetNormal(su2double *val_normal);
  
  /*!
   * \brief Set the value of the volume of the control volume.
   * \param[in] val_volume Volume of the control volume.
   */
  void SetVolume(su2double val_volume);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhosIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhoIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetPIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetTIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetTveIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the velocity index in the primitive variable vector.
   * \param[in] i(rho*u)
   */
  void SetVelIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetHIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetAIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhoCvtrIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhoCvveIndex(unsigned short val_Index);
  
  /*!
   * \brief Sets the value of the derivative of pressure w.r.t. species density.
   * \param[in] iRho_s
   */
  void SetdPdU(su2double *val_dPdU_i, su2double *val_dPdU_j);
  
  /*!
   * \brief Sets the value of the derivative of temperature w.r.t. species density.
   * \param[in] iRho_s
   */
  void SetdTdU(su2double *val_dTdU_i, su2double *val_dTdU_j);
  
  /*!
   * \brief Sets the value of the derivative of vib-el. temperature w.r.t. species density.
   * \param[in] iRho_s
   */
  void SetdTvedU(su2double *val_dTvedU_i, su2double *val_dTvedU_j);
  
  /*!
   * \brief Get the inviscid fluxes.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_pressure - Value of the pressure.
   * \param[in] val_enthalpy - Value of the enthalpy.
   */
  void GetInviscidFlux(su2double val_density, su2double *val_velocity, su2double val_pressure, su2double val_enthalpy);
  
  /*!
   * \brief Get the viscous fluxes.
   * \param[in] val_primvar - Value of the primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_mach_inf - Value of the Mach number at the infinity.
   */
  void GetViscousFlux(su2double *val_primvar, su2double **val_gradprimvar,
                      su2double val_laminar_viscosity, su2double val_eddy_viscosity,
                      su2double val_mach_inf);
  
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
  void GetInviscidArtCompProjFlux(su2double *val_density, su2double *val_velocity,
                                  su2double *val_pressure, su2double *val_betainc2,
                                  su2double *val_normal, su2double *val_Proj_Flux);
  
  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_eddy_conductivity - Eddy Conductivity.
   */
  
  void GetViscousProjFlux(su2double *val_primvar, su2double **val_gradprimvar,
                          su2double val_turb_ke, su2double *val_normal,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity);
  /*!
   * \brief Compute the projection of the viscous fluxes into a direction for general fluid model.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_heat_capacity_cp - Heat Capacity at constant pressure.
   */
  
  void GetViscousProjFlux(su2double *val_primvar, su2double **val_gradprimvar,
                          su2double val_turb_ke, su2double *val_normal,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_heat_capacity_cp);
    
  /*!
   * \brief Compute the projection of the viscous fluxes into a direction for general fluid model.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_reynolds_stress - Reynolds stress.
   * \param[in] val_perturbed_rsm - Perturbed Reynolds Stress Matrix
   */

  void GetViscousProjFlux(su2double *val_primvar, su2double **val_gradprimvar,
                          su2double val_turb_ke, su2double *val_normal,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double **val_reynolds_stress,
                          su2double **val_perturbed_rsm);

  /*
   * \brief Compute the projection of the viscous fluxes into a direction (artificial compresibility method).
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  
  void GetViscousArtCompProjFlux(su2double **val_gradprimvar,
                                 su2double *val_normal,
                                 su2double val_laminar_viscosity,
                                 su2double val_eddy_viscosity);
  
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
   * \brief Compute the projection of the inviscid Jacobian matrices (artificial compresibility).
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidArtCompProjJac(su2double *val_density, su2double *val_velocity,
                                 su2double *val_betainc2, su2double *val_normal,
                                 su2double val_scale,
                                 su2double **val_Proj_Jac_tensor);
  
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
   * \brief TSL-Approximation of Viscous NS Jacobians.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(su2double *val_Mean_PrimVar,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_dist_ij,
                          su2double *val_normal, su2double val_dS,
                          su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);
  
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
                          su2double **val_gradprimvar,
                          su2double *val_Mean_SecVar,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_heat_capacity_cp,
                          su2double val_dist_ij,
                          su2double *val_normal, su2double val_dS,
                          su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);
  
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
   * \brief Compute the projection of the viscous Jacobian matrices.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousArtCompProjJacs(su2double val_laminar_viscosity,
                                 su2double val_eddy_viscosity, su2double val_dist_ij,
                                 su2double *val_normal, su2double val_dS,
                                 su2double **val_Proj_Jac_Tensor_i,
                                 su2double **val_Proj_Jac_Tensor_j);
  
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
                  su2double *val_soundspeed, su2double *val_enthalpy, su2double *val_chi, su2double *val_kappa,
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
        void ComputeResJacobianGiles(CFluidModel *FluidModel, su2double pressure, su2double density, su2double *turboVel, su2double alphaInBC, su2double gammaInBC,  su2double **R_c, su2double **R_c_inv);

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
	 * \param[out]char_jump - pointer to the vector containing the characteristic variable jump.
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
	void GetPrecondJacobian(su2double Beta2, su2double r_hat, su2double s_hat, su2double t_hat, su2double rB2a2, su2double* val_Lambda, su2double* val_normal, su2double** val_absPeJac);
    
	/*!
	 * \brief Computation of the matrix P (artificial compresibility), this matrix diagonalize the conservative Jacobians in
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_betainv2 - Value of the compresibility factor.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPArtCompMatrix(su2double *val_density, su2double *val_velocity,
                         su2double *val_betainv2, su2double *val_normal,
                         su2double **val_p_tensor);
  
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
   * \brief Computation of the matrix P^{-1} (artificial compresibility), this matrix diagonalize the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_betainv2 - Value of the compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPArtCompMatrix_inv(su2double *val_density, su2double *val_velocity,
                             su2double *val_betainv2, su2double *val_normal,
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
  virtual void ComputeResidual(su2double *val_residual, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  virtual void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j);
  
  virtual void ComputeResidual_TransLM(su2double *val_residual,
                                       su2double **val_Jacobian_i,
                                       su2double **val_Jacobian_j, CConfig *config,
                                       su2double &gamma_sep) ;
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual_i,
                               su2double *val_residual_j, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                               su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[out] val_JacobianMeanFlow_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_JacobianMeanFlow_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                               su2double **val_Jacobian_j,
                               su2double **val_JacobianMeanFlow_i,
                               su2double **val_JacobianMeanFlow_j,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_resconv, su2double *val_resvisc,
                               su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                               CConfig *config);
  
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
  virtual void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                               su2double **val_Jacobian_ii,
                               su2double **val_Jacobian_ij,
                               su2double **val_Jacobian_ji,
                               su2double **val_Jacobian_jj, CConfig *config);
  
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
  virtual void ComputeResidual(su2double *val_resconv_i, su2double *val_resvisc_i,
                               su2double *val_resconv_j, su2double *val_resvisc_j,
                               su2double **val_Jacobian_ii,
                               su2double **val_Jacobian_ij,
                               su2double **val_Jacobian_ji,
                               su2double **val_Jacobian_jj, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_stiffmatrix_elem, CConfig *config);
  
  /*!
   * \overload
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] - Matrix for storing the constants to be used in the calculation of the equilibrium extent of reaction Keq.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GetEq_Rxn_Coefficients(su2double **EqnRxnConstants, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual_Axisymmetric(su2double *val_residual, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual_Axisymmetric_ad(su2double *val_residual, su2double *val_residuald, CConfig *config);
  
  /*!
   * \brief Calculation of axisymmetric source term Jacobian
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetJacobian_Axisymmetric(su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculation of the translational-vibrational energy exchange source term
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  virtual void ComputeVibRelaxation(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculation of the chemistry source term
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  virtual void ComputeChemistry(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculates constants used for Keq correlation.
   * \param[out] A - Pointer to coefficient array.
   * \param[in] val_reaction - Reaction number indicator.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GetKeqConstants(su2double *A, unsigned short val_reaction, CConfig *config);
  
  /*!
   * \brief Set intermittency for numerics (used in SA with LM transition model)
   */
  virtual void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  virtual void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  virtual void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  virtual void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  virtual su2double GetProduction(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  virtual su2double GetDestruction(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  virtual su2double GetCrossProduction(void);

  /*!
   * \brief A virtual member.
   */
  virtual su2double GetGammaBC(void);
  
  /*!
   * \overload
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_Jacobian_i,
                               su2double *val_Jacobian_mui,
                               su2double ***val_Jacobian_gradi, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_Jacobian_i,
                               su2double *val_Jacobian_mui,
                               su2double ***val_Jacobian_gradi,
                               su2double **val_Jacobian_j,
                               su2double *val_Jacobian_muj,
                               su2double ***val_Jacobian_gradj, CConfig *config);

  /*!
   * \brief A virtual member to compute the tangent matrix in structural problems
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Tangent_Matrix(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the pressure term in incompressible or nearly-incompressible structural problems
   * \param[in] element_container - Definition of the particular element integrated.
   */
  virtual void Compute_MeanDilatation_Term(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the nodal stress term in non-linear structural problems
   * \param[in] element_container - Definition of the particular element integrated.
   */
  virtual void Compute_NodalStress_Term(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the plane stress term in an element for nonlinear structural problems
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the constitutive matrix in an element for structural problems
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the stress tensor in an element for structural problems
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the element-based Lame parameters and set the local properties
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void SetElement_Properties(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member
   * \param[in] config - Config structure
   */
  virtual void ReadDV(CConfig *config);

  /*!
   * \brief A virtual member to set the value of the design variables
   * \param[in] i_DV - Index of the design variable.
   * \param[in] val_DV - Value of the design variable
   */
  virtual void Set_DV_Val(unsigned short i_DV, su2double val_DV);

  /*!
   * \brief A virtual member to retrieve the value of the design variables
   * \param[in] i_DV - Index of the design variable.
   */
  virtual su2double Get_DV_Val(unsigned short i_DV);

  /*!
   * \brief A virtual member to add the Maxwell stress contribution
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Add_MaxwellStress(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to set element-based electric field modulus
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void SetElectric_Properties(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to set the electric field
   * \param[in] EField_DV - New electric field computed by adjoint methods.
   */
  virtual void Set_ElectricField(unsigned short i_DV, su2double val_EField);

  /*!
   * \brief A virtual member to set the young modulus
   * \param[in] val_Young - Value of the Young Modulus.
   */
  virtual void Set_YoungModulus(unsigned short i_DV, su2double val_Young);

  /*!
   * \brief A virtual member to set the material properties
   * \param[in] iVal - Index of the region of concern
   * \param[in] val_E - Value of the Young Modulus.
   * \param[in] val_Nu - Value of the Poisson's ratio.
   */
  virtual void SetMaterial_Properties(unsigned short iVal, su2double val_E, su2double val_Nu);

  /*!
   * \brief A virtual member to set the material properties
   * \param[in] iVal - Index of the region of concern
   * \param[in] val_Rho - Value of the density (inertial effects).
   * \param[in] val_Rho_DL - Value of the density (dead load effects).
   */
  virtual void SetMaterial_Density(unsigned short iVal, su2double val_Rho, su2double val_Rho_DL);

  /*!
   * \brief A virtual member to compute the mass matrix
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Mass_Matrix(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the residual component due to dead loads
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Dead_Load(CElement *element_container, CConfig *config);

  /*!
   * \brief A virtual member to compute the averaged nodal stresses
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  virtual void Compute_Averaged_NodalStress(CElement *element_container, CConfig *config);

  /*!
   * \brief Computes a basis of orthogonal vectors from a suppled vector
   * \param[in] config - Normal vector
   */
  void CreateBasis(su2double *val_Normal);

  /*!
   * \brief Decomposes the symmetric matrix A_ij, into eigenvectors and eigenvalues
   * \param A_i: symmetric matrix to be decomposed
   * \param Eig_Vec: strores the eigenvectors
   * \param Eig_Val: stores the eigenvalues
   */
  void EigenDecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val);

  /*!
   * \brief Recomposes the eigenvectors and eigenvalues into a matrix
   * \param A_ij: recomposed matrix
   * \param Eig_Vec: eigenvectors
   * \param Eig_Val: eigenvalues
   */
  void EigenRecomposition(su2double **A_ij, su2double **Eig_Vec, su2double *Eig_Val);

  /*!
   * \brief tred2
   * \param V
   * \param d
   * \param e
   */
  void tred2(su2double **V, su2double *d, su2double *e);

  /*!
   * \brief tql2
   * \param V
   * \param d
   * \param e
   */
  void tql2(su2double **V, su2double *d, su2double *e);

};

/*!
 * \class CUpwCUSP_Flow
 * \brief Class for centered scheme - CUSP.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CUpwCUSP_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Flux, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity,  /*!< \brief Mean and projected velocities. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  *ProjFlux, *ProjFlux_i, *ProjFlux_j,  /*!< \brief Projected inviscid flux tensor. */
  cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  LamdaNeg, LamdaPos, Beta, Nu_c, U_i[5], U_j[5], MeanSoundSpeed, Mach,
  **Jacobian;  /*!< \brief Projected grid velocity. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwCUSP_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_residual - Pointer to the residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author A. Bueno, F. Palacios
 * \version 5.0.0 "Raven"
 */
class CUpwRoe_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon, MaxLambda, Delta;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, RoeSoundSpeed2, kappa;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};


/*!
 * \class CUpwGeneralRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations for a general fluid model.
 * \ingroup ConvDiscr
 * \author S.Vitale, G.Gori, M.Pini
 * \version 5.0.0 "Raven"
 */
class CUpwGeneralRoe_Flow : public CNumerics {
private:

  bool implicit, grid_movement;

  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon, MaxLambda, Delta;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,

  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeSoundSpeed2,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho, kappa;
  unsigned short iDim, iVar, jVar, kVar;


  su2double StaticEnthalpy_i, StaticEnergy_i, StaticEnthalpy_j, StaticEnergy_j, Kappa_i, Kappa_j, Chi_i, Chi_j, Velocity2_i, Velocity2_j;
  su2double RoeKappa, RoeChi;

public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwGeneralRoe_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Compute the Average for a general fluid flux between two nodes i and j.
   * Using the approach of Vinokur and Montagne'
   */
  
  void ComputeRoeAverage();
};


/*!
 * \class CUpwMSW_Flow
 * \brief Class for solving a flux-vector splitting method by Steger & Warming, modified version.
 * \ingroup ConvDiscr
 * \author S. Copeland
 * \version 5.0.0 "Raven"
 */
class CUpwMSW_Flow : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double *u_i, *u_j, *ust_i, *ust_j;
  su2double *Fc_i, *Fc_j;
  su2double *Lambda_i, *Lambda_j;
  su2double rhos_i, rhos_j;
  su2double *Ust_i, *Ust_j, *Vst_i, *Vst_j, *Velst_i, *Velst_j;
  su2double **P_Tensor, **invP_Tensor;
  unsigned short nPrimVar, nVar, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwMSW_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
};

/*!
 * \class CUpwTurkel_Flow
 * \brief Class for solving an approximate Riemann solver of Roe with Turkel Preconditioning for the flow equations.
 * \ingroup ConvDiscr
 * \author A. K. Lonkar
 * \version 5.0.0 "Raven"
 */
class CUpwTurkel_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **absPeJac, **invRinvPe, **R_Tensor, **Matrix, **Art_Visc;
  su2double sq_vel, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoePressure, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double Beta, Beta_min, Beta_max;
  su2double r_hat, s_hat, t_hat, rhoB2a2, sqr_one_m_Betasqr_Lam1;
  su2double Beta2, one_m_Betasqr, one_p_Betasqr, sqr_two_Beta_c_Area;
  su2double local_Mach;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwTurkel_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Get the Preconditioning Beta.
   * \return Beta - Value of the low Mach Preconditioner.
   */
  su2double GetPrecond_Beta();
};

/*!
 * \class CUpwArtComp_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the incompressible flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CUpwArtComp_Flow : public CNumerics {
private:
  bool implicit;
  bool gravity;
  su2double Froude;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double Proj_ModJac_Tensor_ij, Pressure_i,
  Pressure_j, MeanDensity, MeanSoundSpeed, MeanPressure, MeanBetaInc2,
  ProjVelocity;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwArtComp_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwRoe_AdjFlow
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CUpwRoe_AdjFlow : public CNumerics {
private:
  su2double *Residual_Roe;
  su2double area, Sx, Sy, Sz, rarea, nx, ny, nz, rho_l, u_l, v_l, w_l, h_l, rho_r,
  u_r, v_r, w_r, h_r, psi1, psi2, psi3, psi4, psi5;
  su2double h, u, v, w, c, psi1_l, psi2_l, psi3_l, psi4_l, psi5_l,
  psi1_r, psi2_r, psi3_r, psi4_r, psi5_r, q_l, q_r, Q_l, Q_r, vn,
  rrho_l, weight, rweight1, cc;
  su2double l1psi, l2psi, absQ, absQp, absQm, q2, alpha, beta_u, beta_v, beta_w, Q, l1l2p, l1l2m, eta;
  su2double RoeDensity, RoeSoundSpeed, *RoeVelocity, *Lambda, *Velocity_i, *Velocity_j, **ProjFlux_i, **ProjFlux_j,
  Proj_ModJac_Tensor_ij, **Proj_ModJac_Tensor, Energy_i, Energy_j, **P_Tensor, **invP_Tensor;
  unsigned short iDim, iVar, jVar, kVar;
  bool implicit, grid_movement;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CUpwRoeArtComp_AdjFlow
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CUpwRoeArtComp_AdjFlow : public CNumerics {
private:
  su2double Area, *Lambda, *Velocity_i, *Velocity_j, **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j,
  Proj_ModJac_Tensor_ij, **Proj_ModJac_Tensor, **P_Tensor, **invP_Tensor, MeanDensity,
  MeanPressure, MeanBetaInc2, ProjVelocity, *MeanVelocity, MeanSoundSpeed;
  unsigned short iDim, iVar, jVar, kVar;
  bool implicit;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoeArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoeArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CUpwAUSM_Flow : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSM_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwHLLC_Flow
 * \brief Class for solving an approximate Riemann HLLC.
 * \ingroup ConvDiscr
 * \author G. Gori, Politecnico di Milano
 * \version 4.0.2 "Cardinal"
 */
class CUpwHLLC_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  unsigned short iDim, jDim, iVar, jVar;
  
  su2double *IntermediateState;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;

  su2double sq_vel_i, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, ProjVelocity_i;
  su2double sq_vel_j, Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, ProjVelocity_j;
  
  su2double sq_velRoe, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeProjVelocity, ProjInterfaceVel;

  su2double sL, sR, sM, pStar, EStar, rhoSL, rhoSR, Rrho, kappa;

  su2double Omega, RHO, OmegaSM;
  su2double *dSm_dU, *dPI_dU, *drhoStar_dU, *dpStar_dU, *dEStar_dU;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwHLLC_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

};

/*!
 * \class CUpwGeneralHLLC_Flow
 * \brief Class for solving an approximate Riemann HLLC.
 * \ingroup ConvDiscr
 * \author G. Gori, Politecnico di Milano
 * \version 4.0.2 "Cardinal"
 */
class CUpwGeneralHLLC_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  unsigned short iDim, jDim, iVar, jVar;
  
  su2double *IntermediateState;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;

  su2double sq_vel_i, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, ProjVelocity_i, StaticEnthalpy_i, StaticEnergy_i;
  su2double sq_vel_j, Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, ProjVelocity_j, StaticEnthalpy_j, StaticEnergy_j;
  
  su2double sq_velRoe, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeProjVelocity, ProjInterfaceVel;
  su2double Kappa_i, Kappa_j, Chi_i, Chi_j, RoeKappa, RoeChi, RoeKappaStaticEnthalpy;

  su2double sL, sR, sM, pStar, EStar, rhoSL, rhoSR, Rrho, kappa;

  su2double Omega, RHO, OmegaSM;
  su2double *dSm_dU, *dPI_dU, *drhoStar_dU, *dpStar_dU, *dEStar_dU;

  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwGeneralHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwGeneralHLLC_Flow(void);
  
  /*!

   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
   void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

   /*!
   * \brief Compute the Average quantities for a general fluid flux between two nodes i and j.
   * Using the approach of Vinokur and Montagne'
   */
   void VinokurMontagne();
};

/*!
 * \class CUpwLin_TransLM
 * \brief Class for performing a linear upwind solver for the Spalart-Allmaras turbulence model equations with transition
 * \ingroup ConvDiscr
 * \author A. Aranake
 * \version 5.0.0 "Raven"
 */
class CUpwLin_TransLM : public CNumerics {
private:
  su2double *Velocity_i;
  su2double *Velocity_j;
  bool implicit, grid_movement, incompressible;
  su2double Density_i, Density_j, q_ij, a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLin_TransLM(void);
  
  /*!
   * \brief Compute the upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwLin_AdjTurb
 * \brief Class for performing a linear upwind solver for the adjoint turbulence equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CUpwLin_AdjTurb : public CNumerics {
private:
  su2double *Velocity_i;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLin_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLin_AdjTurb(void);
  
  /*!
   * \brief Compute the adjoint upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

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
 * \version 5.0.0 "Raven"
 */
class CUpwScalar : public CNumerics {
private:

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(su2double *val_residual,
                                  su2double **Jacobian_i,
                                  su2double **Jacobian_j,
                                  CConfig *config) = 0;

protected:
  su2double *Velocity_i, *Velocity_j; /*!< \brief Velocity, minus any grid movement. */
  su2double Density_i, Density_j;
  bool implicit, grid_movement, incompressible;
  su2double q_ij, /*!< \brief Projected velocity at the face. */
            a0,   /*!< \brief The maximum of the face-normal velocity and 0 */
            a1;   /*!< \brief The minimum of the face-normal velocity and 0 */
  unsigned short iDim;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwScalar(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwScalar(void);

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaras turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CUpwSca_TurbSA : public CUpwScalar {
private:

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn();

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TurbSA(void);
};

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 * \version 5.0.0 "Raven"
 */
class CUpwSca_TurbSST : public CUpwScalar {
private:

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn();

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TurbSST(void);
};

/*!
 * \class CUpwSca_TransLM
 * \brief Class for doing a scalar upwind solver for the Spalart-Allmaras turbulence model equations with transition.
 * \ingroup ConvDiscr
 * \author A. Aranake.
 * \version 5.0.0 "Raven"
 */
class CUpwSca_TransLM : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  bool implicit, grid_movement;
  su2double q_ij, a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TransLM(void);
  
  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_AdjTurb
 * \brief Class for doing a scalar upwind solver for the adjoint turbulence equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CUpwSca_AdjTurb : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_AdjTurb(void);
  
  /*!
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};


/*!
 * \class CCentJST_Flow
 * \brief Class for centered shceme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentJST_KE_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  ProjGridVel;  /*!< \brief Projected grid velocity. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_KE_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentJST_Flow
 * \brief Class for centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentJST_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, Epsilon_4, cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  ProjGridVel;  /*!< \brief Projected grid velocity. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentJSTArtComp_Flow
 * \brief Class for centered scheme - JST (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentJSTArtComp_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanBetaInc2, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, Epsilon_4, cte_0, cte_1; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement, /*!< \brief Modification for grid movement. */
  gravity; /*!< \brief computation with gravity force. */
  su2double Froude; /*!< \brief Froude number. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJSTArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJSTArtComp_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentJST_AdjFlow
 * \brief Class for and adjoint centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentJST_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi;
  unsigned short iDim, jDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2;
  su2double MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_4, Param_Kappa_2, Local_Lambda_i, Local_Lambda_j, MeanLambda;
  su2double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_2;
  bool implicit, grid_movement;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a JST method.
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
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentJSTArtComp_AdjFlow
 * \brief Class for and adjoint centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentJSTArtComp_AdjFlow : public CNumerics {
private:
  su2double sc2, *Diff_Psi, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j;
  unsigned short iDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j;
  su2double Param_p, Param_Kappa_4, Param_Kappa_2, Local_Lambda_i, Local_Lambda_j, MeanLambda;
  su2double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_2;
  bool implicit;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJSTArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJSTArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a JST method.
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
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentLax_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentLax_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, /*!< \brief Difference of conservative variables. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
  Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  su2double ProjGridVel;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentLaxArtComp_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentLaxArtComp_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, /*!< \brief Difference of conservative variables. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanBetaInc2, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
  Epsilon_0; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement, /*!< \brief Modification for grid movement. */
  gravity; /*!< \brief Modification for for gravity force. */
  su2double Froude;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLaxArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLaxArtComp_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentLax_AdjFlow
 * \brief Class for computing the Lax-Friedrich adjoint centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentLax_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi;
  unsigned short iDim, jDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2,
  MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_0, Local_Lambda_i, Local_Lambda_j, MeanLambda,
  Phi_i, Phi_j, sc2, StretchingFactor, Epsilon_0;
  bool implicit, grid_movement;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a Lax method.
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
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentLaxArtComp_AdjFlow
 * \brief Class for computing the Lax-Friedrich adjoint centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CCentLaxArtComp_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi, **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j;
  unsigned short iDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, Param_p, Param_Kappa_0,
  Local_Lambda_i, Local_Lambda_j, MeanLambda,
  Phi_i, Phi_j, sc2, StretchingFactor, Epsilon_0;
  bool implicit;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLaxArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLaxArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a Lax method.
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
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;     /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,           /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,           /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,             /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij;            /*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
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
 * \version 3.2.1 "eagle"
 */

class CGeneralAvgGrad_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;     /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,           /*!< \brief Mean primitive variables. */
  *Mean_SecVar,                   /*!< \brief Mean secondary variables. */
  *PrimVar_i, *PrimVar_j,           /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,             /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,             /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp,                               /*!< \brief Mean value of the Cp. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij;            /*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
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
 * \class CAvgGradArtComp_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGradArtComp_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;  /*!< \brief Iterators in dimension an variable. */
  su2double **Mean_GradPrimVar,          /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
  dist_ij;              /*!< \brief Length of the edge and face. */
  bool implicit;        /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradArtComp_Flow(void);
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
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_Scalar : public CNumerics {
 private:

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(su2double *val_residual,
                                  su2double **Jacobian_i,
                                  su2double **Jacobian_j,
                                  CConfig *config) = 0;

 protected:
  bool implicit, incompressible;
  bool correct_gradient;
  unsigned short iVar, iDim;
  su2double **Mean_GradTurbVar;               /*!< \brief Average of gradients at cell face */
  su2double *Edge_Vector,                     /*!< \brief Vector from node i to node j. */
            *Proj_Mean_GradTurbVar_Normal,    /*!< \brief Mean_gradTurbVar DOT normal */
            *Proj_Mean_GradTurbVar_Edge,      /*!< \brief Mean_gradTurbVar DOT Edge_Vector */
            *Proj_Mean_GradTurbVar;           /*!< \brief Mean_gradTurbVar DOT normal, corrected if required*/
  su2double  dist_ij_2,                       /*!< \brief |Edge_Vector|^2 */
             proj_vector_ij;                  /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Scalar(unsigned short val_nDim, unsigned short val_nVar,
                    bool correct_gradient, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Scalar(void);

  /*!
   * \brief Compute the viscous residual using an average of gradients without correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                       su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_TurbSA : public CAvgGrad_Scalar {
private:

  const su2double sigma;
  su2double nu_i, nu_j, nu_e;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSA(void);
};

/*!
 * \class CAvgGrad_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_TurbSA_Neg : public CAvgGrad_Scalar {
private:

  const su2double sigma;
  const su2double cn1;
  su2double fn, Xi;
  su2double nu_i, nu_j, nu_ij, nu_tilde_ij, nu_e;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                      bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSA_Neg(void);
};

/*!
 * \class CAvgGrad_TransLM
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_TransLM : public CNumerics {
private:
  su2double **Mean_GradTransVar;
  su2double *Proj_Mean_GradTransVar_Kappa, *Proj_Mean_GradTransVar_Edge;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  //su2double dist_ij_2;
  //su2double proj_vector_ij;
  //unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TransLM(void);
  
  /*!
   * \brief Compute the viscous turbulence terms residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;  /*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;  /*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double *Mean_GradPsiE;  /*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;  /*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;  /*!< \brief Vector going from node i to node j. */
  bool implicit;      /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjFlow(void);
  
  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradArtComp_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGradArtComp_AdjFlow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;  /*!< \brief Iterators in dimension an variable. */
  su2double **Mean_GradPsiVar,          /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
  dist_ij;              /*!< \brief Length of the edge and face. */
  bool implicit;        /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradArtComp_AdjFlow(void);
  
  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;    /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,          /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,        /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                  /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,         /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2;           /*!< \brief Length of the edge and face. */
  bool implicit;      /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

  /*!
   * \brief Initialize the Reynolds Stress Matrix
   * \param turb_ke turbulent kinetic energy of node
   */
  void SetReynoldsStressMatrix(su2double turb_ke);

  /*!
   * \brief Perturb the Reynolds stress tensor based on parameters
   * \param turb_ke: turbulent kinetic energy of the noce
   * \param Eig_Val_Comp: Defines type of eigenspace perturbation
   * \param beta_delta: Defines the amount of eigenvalue perturbation
   */
  void SetPerturbedRSM(su2double turb_ke, CConfig *config);

  /*!
   * \brief Get the mean rate of strain matrix based on velocity gradients
   * \param S_ij
   */
  void GetMeanRateOfStrainMatrix(su2double **S_ij);

  /*!
   * \brief Setting the UQ framework usage
   * \param val_using_uq
   */
  void SetUsing_uq(bool val_using_uq);
};


/*!
 * \class CGeneralAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author M. Pini, S. Vitale
 * \version 3.2.1 "eagle"
 */
class CGeneralAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;    /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,          /*!< \brief Mean primitive variables. */
  *Mean_SecVar,                  /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,            /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                  /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,   /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp,                     /*!< \brief Mean value of the specific heat. */
  Mean_turb_ke,         /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2;           /*!< \brief Length of the edge and face. */
  bool implicit;      /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGradCorrected_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrectedArtComp_Flow
 * \brief Class for computing viscous term using an average of gradients with correction (artificial compresibility).
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGradCorrectedArtComp_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;  /*!< \brief Iterators in dimension an variable. */
  su2double *PrimVar_i, *PrimVar_j,      /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,      /*!< \brief Mean value of the viscosity. */
  dist_ij_2;          /*!< \brief Length of the edge and face. */
  bool implicit;      /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrectedArtComp_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_TransLM
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CAvgGradCorrected_TransLM : public CNumerics {
private:
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TransLM(void);
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_TurbSST : public CAvgGrad_Scalar {
private:
  su2double sigma_k1, /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_k2,
  sigma_om1,
  sigma_om2;

  su2double diff_kine,  /*!< \brief Diffusivity for viscous terms of tke eq */
            diff_omega; /*!< \brief Diffusivity for viscous terms of omega eq */

  su2double F1_i, F1_j; /*!< \brief Menter's first blending function */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                   su2double* constants, bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSST(void);

  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) {
    F1_i = val_F1_i; F1_j = val_F1_j;
  }

};

/*!
 * \class CAvgGradCorrected_AdjFlow
 * \brief Class for computing the adjoint viscous terms, including correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CAvgGradCorrected_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;  /*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;  /*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double **Mean_GradPsiVar;  /*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;  /*!< \brief Vector going from node i to node j. */
  su2double *Proj_Mean_GradPsiVar_Edge;  /*!< \brief Projection of Mean_GradPsiVar onto Edge_Vector. */
  su2double *Mean_GradPsiE;  /*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;  /*!< \brief Counter for dimensions of the problem. */
  bool implicit;          /*!< \brief Boolean controlling Jacobian calculations. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
   * \param[out] val_residual_i - Pointer to the viscous residual at point i.
   * \param[out] val_residual_j - Pointer to the viscous residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradCorrectedArtComp_AdjFlow
 * \brief Class for computing the adjoint viscous terms, including correction.
 * \ingroup ViscDiscr
 * \author F.Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGradCorrectedArtComp_AdjFlow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;  /*!< \brief Iterators in dimension an variable. */
  su2double *PsiVar_i, *PsiVar_j,      /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                /*!< \brief Vector form point i to point j. */
  **Mean_GradPsiVar, *Proj_Mean_GradPsiVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,      /*!< \brief Mean value of the viscosity. */
  dist_ij_2;          /*!< \brief Length of the edge and face. */
  bool implicit;      /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrectedArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrectedArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
   * \param[out] val_residual_i - Pointer to the viscous residual at point i.
   * \param[out] val_residual_j - Pointer to the viscous residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                        su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CAvgGradCorrected_AdjTurb : public CNumerics {
private:
  su2double **Mean_GradTurbPsi;
  su2double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
  su2double *Edge_Vector;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_AdjTurb(void);
  
  /*!
   * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
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
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGrad_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CAvgGrad_AdjTurb : public CNumerics {
private:
  su2double **Mean_GradTurbPsi;
  su2double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
  su2double *Edge_Vector;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjTurb(void);
  
  /*!
   * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
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
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CGalerkin_Flow
 * \brief Class for computing the stiffness matrix of the Galerkin method.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CGalerkin_Flow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGalerkin_Flow(void);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double **val_stiffmatrix_elem, CConfig *config);
};

/*!
 * \class CFEM_Elasticity
 * \brief Generic class for computing the tangent matrix and the residual for structural problems
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_Elasticity : public CNumerics {

protected:

  su2double E;              /*!< \brief Aux. variable, Young's modulus of elasticity. */
  su2double Nu;             /*!< \brief Aux. variable, Poisson's ratio. */
  su2double Rho_s;          /*!< \brief Aux. variable, Structural density. */
  su2double Rho_s_DL;       /*!< \brief Aux. variable, Structural density (for dead loads). */
  su2double Mu;             /*!< \brief Aux. variable, Lame's coeficient. */
  su2double Lambda;         /*!< \brief Aux. variable, Lame's coeficient. */
  su2double Kappa;          /*!< \brief Aux. variable, Compressibility constant. */

  su2double *E_i;           /*!< \brief Young's modulus of elasticity. */
  su2double *Nu_i;          /*!< \brief Poisson's ratio. */
  su2double *Rho_s_i;       /*!< \brief Structural density. */
  su2double *Rho_s_DL_i;    /*!< \brief Structural density (for dead loads). */

  bool plane_stress;        /*!< \brief Checks if we are solving a plane stress case */

  su2double **Ba_Mat,          /*!< \brief Matrix B for node a - Auxiliary. */
  **Bb_Mat;                    /*!< \brief Matrix B for node b - Auxiliary. */
  su2double *Ni_Vec;           /*!< \brief Vector of shape functions - Auxiliary. */
  su2double **D_Mat;           /*!< \brief Constitutive matrix - Auxiliary. */
  su2double **KAux_ab;         /*!< \brief Node ab stiffness matrix - Auxiliary. */
  su2double **GradNi_Ref_Mat;  /*!< \brief Gradients of Ni - Auxiliary. */
  su2double **GradNi_Curr_Mat; /*!< \brief Gradients of Ni - Auxiliary. */

  su2double *FAux_Dead_Load;    /*!< \brief Auxiliar vector for the dead loads */

  su2double *DV_Val;          /*!< \brief For optimization cases, value of the design variables. */
  unsigned short n_DV;          /*!< \brief For optimization cases, number of design variables. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_Elasticity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFEM_Elasticity(void);

  void SetMaterial_Properties(unsigned short iVal, su2double val_E, su2double val_Nu);

  void SetMaterial_Density(unsigned short iVal, su2double val_Rho, su2double val_Rho_DL);

  void Compute_Mass_Matrix(CElement *element_container, CConfig *config);

  void Compute_Dead_Load(CElement *element_container, CConfig *config);

  void Set_YoungModulus(unsigned short i_DV, su2double val_Young);

  void SetElement_Properties(CElement *element_container, CConfig *config);

  void ReadDV(CConfig *config);

  void Set_DV_Val(unsigned short i_DV, su2double val_DV);

  su2double Get_DV_Val(unsigned short i_DV);

  virtual void Compute_Tangent_Matrix(CElement *element_container, CConfig *config);

  virtual void Compute_MeanDilatation_Term(CElement *element_container, CConfig *config);

  virtual void Compute_NodalStress_Term(CElement *element_container, CConfig *config);

  virtual void Compute_Averaged_NodalStress(CElement *element_container, CConfig *config);

  virtual void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  virtual void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);
  
  virtual void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

	virtual void Add_MaxwellStress(CElement *element_container, CConfig *config);

  virtual void SetElectric_Properties(CElement *element_container, CConfig *config);

  virtual void Set_ElectricField(unsigned short i_DV, su2double val_EField);

};

/*!
 * \class CFEM_LinearElasticity
 * \brief Class for computing the stiffness matrix of a linear, elastic problem.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_LinearElasticity : public CFEM_Elasticity {

  su2double **nodalDisplacement;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_LinearElasticity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_LinearElasticity(void);

  void Compute_Tangent_Matrix(CElement *element_container, CConfig *config);

  void Compute_Constitutive_Matrix(void);
  using CNumerics::Compute_Constitutive_Matrix;

  void Compute_Averaged_NodalStress(CElement *element_container, CConfig *config);

};

/*!
 * \class CFEM_LinearElasticity
 * \brief Class for computing the stiffness matrix of a nonlinear, elastic problem.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_NonlinearElasticity : public CFEM_Elasticity {

protected:

  su2double **F_Mat;             /*!< \brief Deformation gradient. */
  su2double **b_Mat;             /*!< \brief Left Cauchy-Green Tensor. */
  su2double **currentCoord;      /*!< \brief Current coordinates. */
  su2double **Stress_Tensor;     /*!< \brief Cauchy stress tensor */

  su2double **FmT_Mat;           /*!< \brief Deformation gradient inverse and transpose. */

  su2double **KAux_P_ab;         /*!< \brief Auxiliar matrix for the pressure term */
  su2double *KAux_t_a;           /*!< \brief Auxiliar matrix for the pressure term */

  su2double J_F;                 /*!< \brief Jacobian of the transformation (determinant of F) */

  su2double f33;                 /*!< \brief Plane stress term for non-linear 2D plane stress analysis */

  bool nearly_incompressible;    /*!< \brief Boolean to consider nearly_incompressible effects */
  bool incompressible;           /*!< \brief Boolean to consider Hu-Washizu incompressible effects */

  su2double **F_Mat_Iso;         /*!< \brief Isocoric component of the deformation gradient. */
  su2double **b_Mat_Iso;         /*!< \brief Isocoric component of the left Cauchy-Green tensor. */

  su2double C10, D1;             /*!< \brief C10 = Mu/2. D1 = Kappa/2. */
  su2double J_F_Iso;             /*!< \brief J_F_Iso: det(F)^-1/3. */

  su2double ****cijkl;           /*!< \brief Constitutive tensor i,j,k,l (defined only for incompressibility - near inc.). */

  bool maxwell_stress;           /*!< \brief Consider the effects of the dielectric loads */

  su2double *EField_Ref_Unit,    /*!< \brief Electric Field, unitary, in the reference configuration. */
  *EField_Ref_Mod;               /*!< \brief Electric Field, modulus, in the reference configuration. */
  su2double *EField_Curr_Unit;   /*!< \brief Auxiliary vector for the unitary Electric Field in the current configuration. */
  unsigned short nElectric_Field,
  nDim_Electric_Field;

  su2double *ke_DE_i;           /*!< \brief Electric Constant for Dielectric Elastomers. */

  su2double ke_DE;              /*!< \brief Electric Constant for Dielectric Elastomers. */
  su2double EFieldMod_Ref;      /*!< \brief Modulus of the electric field in the reference configuration. */


public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_NonlinearElasticity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFEM_NonlinearElasticity(void);

  void Compute_Tangent_Matrix(CElement *element_container, CConfig *config);

  void Compute_MeanDilatation_Term(CElement *element_container, CConfig *config);

  void Compute_NodalStress_Term(CElement *element_container, CConfig *config);

  void Compute_Averaged_NodalStress(CElement *element_container, CConfig *config);

  void Add_MaxwellStress(CElement *element_container, CConfig *config);

  void SetElectric_Properties(CElement *element_container, CConfig *config);

  void Compute_FmT_Mat(void);

  void Compute_Isochoric_F_b(void);

  void Assign_cijkl_D_Mat(void);

  void Set_ElectricField(unsigned short i_DV, su2double val_EField);

  void Set_YoungModulus(unsigned short i_DV, su2double val_Young);

  void SetMaterial_Properties(unsigned short iVal, su2double val_E, su2double val_Nu);

  void SetMaterial_Density(unsigned short iVal, su2double val_Rho, su2double val_Rho_DL);

  su2double deltaij(unsigned short iVar, unsigned short jVar);

  virtual void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  virtual void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);

  virtual void Compute_Stress_Tensor(CElement *element_container, CConfig *config);


};

/*!
 * \class CFEM_NeoHookean_Comp
 * \brief Class for computing the constitutive and stress tensors for a neo-Hookean material model, compressible.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_NeoHookean_Comp : public CFEM_NonlinearElasticity {

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_NeoHookean_Comp(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_NeoHookean_Comp(void);

  void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);
  using CNumerics::Compute_Constitutive_Matrix;

  void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

};

/*!
 * \class CFEM_NeoHookean_Incomp
 * \brief Class for computing the constitutive and stress tensors for a neo-Hookean material model, incompressible.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_NeoHookean_Incomp : public CFEM_NonlinearElasticity {

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_NeoHookean_Incomp(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_NeoHookean_Incomp(void);

  void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);
  using CNumerics::Compute_Constitutive_Matrix;

  void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

};

/*!
 * \class CFEM_IdealDE
 * \brief Class for computing the constitutive and stress tensors for a nearly-incompressible ideal DE.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_IdealDE : public CFEM_NonlinearElasticity {

	su2double trbbar, Eg, Eg23, Ek, Pr;	/*!< \brief Variables of the model calculation. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_IdealDE(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_IdealDE(void);

  void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);

  void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

};

/*!
 * \class CFEM_NeoHookean_Comp
 * \brief Class for computing the constitutive and stress tensors for a Knowles stored-energy function, nearly incompressible.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_Knowles_NearInc : public CFEM_NonlinearElasticity {

	su2double trbbar, term1, term2, Ek, Pr;	/*!< \brief Variables of the model calculation. */
	su2double Bk, Nk;						/*!< \brief Parameters b and n of the model. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_Knowles_NearInc(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_Knowles_NearInc(void);

  void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);
  using CNumerics::Compute_Constitutive_Matrix;

	void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

};

/*!
 * \class CFEM_DielectricElastomer
 * \brief Class for computing the constitutive and stress tensors for a dielectric elastomer.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CFEM_DielectricElastomer : public CFEM_NonlinearElasticity {

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_DielectricElastomer(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DielectricElastomer(void);

  void Compute_Plane_Stress_Term(CElement *element_container, CConfig *config);

  void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config);
  using CNumerics::Compute_Constitutive_Matrix;

  void Compute_Stress_Tensor(CElement *element_container, CConfig *config);

};


/*!
 * \class CSourceNothing
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CSourceNothing : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceNothing(void);
};

/*!
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CSourcePieceWise_TurbSA : public CNumerics {
private:
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double ct3;
  su2double ct4;
  su2double cw3_6;
  su2double cb2_sigma;
  su2double sigma;
  su2double cb2;
  su2double cw1;
  unsigned short iDim;
  su2double nu, Ji, fv1, fv2, ft2, Omega, S, Shat, inv_Shat, dist_i_2, Ji_2, Ji_3, inv_k2_d2;
  su2double r, g, g_6, glim, fw;
  su2double norm2_Grad;
  su2double dfv1, dfv2, dShat;
  su2double dr, dg, dfw;
  bool incompressible;
  bool rotating_frame;
  bool transition;
  su2double gamma_BC;
  su2double intermittency;
  su2double Production, Destruction, CrossProduction;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbSA(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief ______________.
   */
  su2double GetProduction(void);

  /*!
   * \brief  Get the intermittency for the BC trans. model.
   * \return Value of the intermittency.
   */
  su2double GetGammaBC(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetDestruction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetCrossProduction(void);
};

/*!
 * \class CSourcePieceWise_TurbSA_Neg
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CSourcePieceWise_TurbSA_Neg : public CNumerics {
private:
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double ct3;
  su2double ct4;
  su2double cw3_6;
  su2double cb2_sigma;
  su2double sigma;
  su2double cb2;
  su2double cw1;
  unsigned short iDim;
  su2double nu, Ji, fv1, fv2, ft2, Omega, S, Shat, inv_Shat, dist_i_2, Ji_2, Ji_3, inv_k2_d2;
  su2double r, g, g_6, glim, fw;
  su2double norm2_Grad;
  su2double dfv1, dfv2, dShat;
  su2double dr, dg, dfw;
  bool incompressible;
  bool rotating_frame;
  su2double intermittency;
  su2double Production, Destruction, CrossProduction;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbSA_Neg(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief ______________.
   */
  su2double GetProduction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetDestruction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetCrossProduction(void);
};

/*!
 * \class CSourcePieceWise_TransLM
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CSourcePieceWise_TransLM : public CNumerics {
private:
  
  /*-- SA model constants --*/
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double cw3_6;
  su2double sigma;
  su2double cb2;
  su2double cw1;
  
  /*-- gamma-theta model constants --*/
  su2double c_e1;
  su2double c_a1;
  su2double c_e2;
  su2double c_a2;
  su2double sigmaf;
  su2double s1;
  su2double c_theta;
  su2double sigmat;
  
  /*-- Correlation constants --*/
  su2double flen_global;
  su2double alpha_global;
  su2double Vorticity;

  bool implicit;
  
public:
  bool debugme; // For debugging only, remove this. -AA
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TransLM(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual_TransLM(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config, su2double &gamma_sep);
  
  void CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(su2double *TransVar_i, su2double *TransVar_id, su2double *val_residual, su2double *val_residuald, CConfig *config);
};

/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 * \version 5.0.0 "Raven"
 */
class CSourcePieceWise_TurbSST : public CNumerics {
private:
  su2double F1_i,
  F1_j,
  F2_i,
  F2_j;
  
  su2double alfa_1,
  alfa_2,
  beta_1,
  beta_2,
  sigma_omega_1,
  sigma_omega_2,
  beta_star,
  a1;
  
  su2double CDkw_i, CDkw_j;
  
  bool incompressible;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double* constants, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbSST(void);
  
  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first blending function at point i.
   * \param[in] val_F1_j - Value of the first blending function at point j.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j);
  
  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F2_i - Value of the second blending function at point i.
   * \param[in] val_F2_j - Value of the second blending function at point j.
   */
  void SetF2blending(su2double val_F2_i, su2double val_F2_j);
  
  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
   */
  virtual void SetCrossDiff(su2double val_CDkw_i, su2double val_CDkw_j);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Initialize the Reynolds Stress Matrix
   * \param turb_ke turbulent kinetic energy of node
   */
  void SetReynoldsStressMatrix(su2double turb_ke);

  /*!
   * \brief Perturb the Reynolds stress tensor based on parameters
   * \param turb_ke: turbulent kinetic energy of the noce
   * \param config: config file
   */
  void SetPerturbedRSM(su2double turb_ke, CConfig *config);
  /*!
     * \brief A virtual member. Get strain magnitude based on perturbed reynolds stress matrix
     * \param turb_ke: turbulent kinetic energy of the node
     */
  void SetPerturbedStrainMag(su2double turb_ke);

  /*!
   * \brief Get the mean rate of strain matrix based on velocity gradients
   * \param S_ij
   */
  void GetMeanRateOfStrainMatrix(su2double **S_ij);

  /*!
   * \brief Setting the UQ framework usage
   * \param val_using_uq
   */
  void SetUsing_uq(bool val_using_uq);
};

/*!
 * \class CSourceGravity
 * \brief Class for the source term integration of the gravity force.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CSourceGravity : public CNumerics {
  su2double Froude;
  bool compressible, incompressible;
  
public:
  
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceGravity(void);
  
  /*!
   * \brief Source term integration for the poissonal potential.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};

/*!
 * \class CSourceBodyForce
 * \brief Class for the source term integration of a body force.
 * \ingroup SourceDiscr
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CSourceBodyForce : public CNumerics {
  su2double *Body_Force_Vector;
  bool compressible;

public:

  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBodyForce(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceBodyForce(void);

  /*!
   * \brief Source term integration for a body force.
   * \param[out] val_residual - Pointer to the residual vector.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);

};

/*!
 * \class CSourceViscous_AdjFlow
 * \brief Class for source term integration in adjoint problem.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CSourceViscous_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *GradDensity, *GradInvDensity, *dPoDensity2, *alpha, *beta, *Sigma_5_vec;
  su2double **GradVel_o_Rho, **sigma, **Sigma_phi, **Sigma_5_Tensor, **Sigma;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceViscous_AdjFlow(void);
  
  /*!
   * \brief Source term integration of the flow adjoint equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, CConfig *config);
  
};

/*!
 * \class CSourcePieceWise_AdjTurb
 * \brief Class for source term integration of the adjoint turbulent equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CSourcePieceWise_AdjTurb : public CNumerics {
private:
  su2double **tau, *Velocity;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_AdjTurb(void);
  
  /*!
   * \brief Source term integration of the adjoint turbulence equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

class CSourceConservative_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *Residual_i, *Residual_j, *Mean_Residual;
  su2double **Mean_PrimVar_Grad;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceConservative_AdjFlow(void);
  
  /*!
   * \brief Source term integration using a conservative scheme.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};

/*!
 * \class CSourceConservative_AdjTurb
 * \brief Class for source term integration in adjoint turbulent problem using a conservative scheme.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 5.0.0 "Raven"
 */
class CSourceConservative_AdjTurb : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceConservative_AdjTurb(void);
  
  /*!
   * \brief Source term integration using a conservative scheme.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CSourceRotatingFrame_Flow
 * \brief Class for a rotating frame source term.
 * \ingroup SourceDiscr
 * \author F. Palacios, T. Economon.
 * \version 5.0.0 "Raven"
 */
class CSourceRotatingFrame_Flow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceRotatingFrame_Flow(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};

/*!
 * \class CSourceRotatingFrame_AdjFlow
 * \brief Source term class for rotating frame adjoint.
 * \ingroup SourceDiscr
 * \author T. Economon.
 * \version 5.0.0 "Raven"
 */
class CSourceRotatingFrame_AdjFlow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceRotatingFrame_AdjFlow(void);
  
  /*!
   * \brief Residual of the adjoint rotating frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};

/*!
 * \class CSourceAxisymmetric_Flow
 * \brief Class for source term for solving axisymmetric problems.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CSourceAxisymmetric_Flow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_Flow(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);
  
};

/*!
 * \class CSourceAxisymmetric_AdjFlow
 * \brief Class for source term for solving axisymmetric problems.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CSourceAxisymmetric_AdjFlow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_AdjFlow(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);
};

/*!
 * \class CSourceWindGust
 * \brief Class for a source term due to a wind gust.
 * \ingroup SourceDiscr
 * \author S. Padrón
 * \version 5.0.0 "Raven"
 */
class CSourceWindGust : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceWindGust(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceWindGust(void);
  
  /*!
   * \brief Residual of the wind gust source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};

/*!
 * \class CSource_Template
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author A. Lonkar.
 * \version 5.0.0 "Raven"
 */
class CSource_Template : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config -  Name of the input config file
   *
   */
  CSource_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSource_Template(void);
};

/*!
 * \class CConvectiveTemplate
 * \brief Class for setting up new method for spatial discretization of convective terms in flow Equations
 * \ingroup ConvDiscr
 * \author A. Lonkar
 * \version 5.0.0 "Raven"
 */
class CConvective_Template : public CNumerics {
private:
  
  /* define private variables here */
  bool implicit;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CConvective_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CConvective_Template(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CViscous_Template
 * \brief Class for computing viscous term using average of gradients.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CViscous_Template : public CNumerics {
private:
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CViscous_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CViscous_Template(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

#include "numerics_structure.inl"
