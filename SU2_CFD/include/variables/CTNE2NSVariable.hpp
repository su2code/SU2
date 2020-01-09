/*!
 * \file CNSVariable.hpp
 * \brief Class for defining the variables of the compressible Navier-Stokes solver.
 * \author F. Palacios, T. Economon, W. Maier, S.R. Copeland
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "CTNE2EulerVariable.hpp"

/*!
 * \class CTNE2NSVariable
 * \brief Main class for defining the variables of the TNE2 Navier-Stokes' solver.
 * \ingroup Navier_Stokes_Equations
 * \author S. R. Copeland, F. Palacios, W. Maier.
 * \version 6.2.0
 */
class CTNE2NSVariable final : public CTNE2EulerVariable {
private:
  VectorType Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
  VectorType Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
  VectorType Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
  VectorType Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
  MatrixType DiffusionCoeff;    /*!< \brief Diffusion coefficient of the mixture. */
  VectorOfMatrix Dij;            /*!< \brief Binary diffusion coefficients. */
  VectorType LaminarViscosity;  /*!< \brief Viscosity of the fluid. */
  VectorType ThermalCond;       /*!< \brief T-R thermal conductivity of the gas mixture. */
  VectorType ThermalCond_ve;    /*!< \brief V-E thermal conductivity of the gas mixture. */

  su2double inv_TimeScale;      /*!< \brief Inverse of the reference time scale. */

  MatrixType Vorticity;         /*!< \brief Vorticity of the fluid. */
  VectorType StrainMag;         /*!< \brief Magnitude of rate of strain tensor. */
  VectorType Tau_Wall;          /*!< \brief Magnitude of the wall shear stress from a wall function. */
  VectorType DES_LengthScale;   /*!< \brief DES Length Scale. */
  VectorType Roe_Dissipation;   /*!< \brief Roe low dissipation coefficient. */
  VectorType Vortex_Tilting;    /*!< \brief Value of the vortex tilting variable for DES length scale computation. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nVarPrim - Number of primitive variables.
   * \param[in] val_nVarPrimGrad - Number of primitive gradient variables.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2NSVariable(unsigned long val_nDim, unsigned long val_nVar,
                  unsigned long val_nPrimVar, unsigned long val_nPrimVarGrad,
                  unsigned long npoint, CConfig *config);

  /*!
   * \brief Constructor of the class.
   * \param[in] val_density - Value of the flow density (initialization value).
   * \param[in] val_massfrac - Value of the flow mass fraction (initialization value).
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_temperature - Value of the flow temperature (initialization value).
   * \param[in] val_temperature_ve - Value of the flow temperature_ve (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nPrimVar - Number of primitive variables.
   * \param[in] val_nPrimVargrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2NSVariable(su2double val_density, su2double *val_massfrac, su2double *val_velocity,
                  su2double val_temperature, su2double val_temperature_ve, unsigned long npoint,
                  unsigned long val_nDim, unsigned long val_nVar, unsigned long val_nPrimVar,
                  unsigned long val_nPrimVarGrad, CConfig *config);

  /*!
   * \brief Constructor of the class.
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nPrimVar - Number of primitive variables.
   * \param[in] val_nPrimgVarGrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2NSVariable(su2double *val_solution, unsigned long val_nDim, unsigned long val_nVar,
                  unsigned long val_nPrimVar, unsigned long val_nPrimVarGrad, unsigned long npoint,
                  CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTNE2NSVariable() = default;

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar_Compressible(unsigned long iPoint, CConfig *config) override;

  /*!
   * \brief Set the laminar viscosity.
   */
  void SetDiffusionCoeff_GuptaYos(CConfig *config);

  /*!
   * \brief Set the laminar viscosity.
   */
  void SetLaminarViscosity_GuptaYos(CConfig *config);

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  void SetThermalConductivity_GuptaYos(CConfig *config);

  /*!
   * \brief Set the vorticity value.
   */
  bool SetVorticity(void);

  /*!
   * \brief Set the transport coefficients for the Wilke/Blottner/Eucken model
   */
  void SetTransportCoefficients_WBE(CConfig *config);

  /*!
   * \brief Get the species diffusion coefficient.
   * \return Value of the species diffusion coefficient.
   */
  inline su2double* GetDiffusionCoeff(unsigned long iPoint) override { return DiffusionCoeff[iPoint]; }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(unsigned long iPoint) const override { return LaminarViscosity(iPoint); }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(unsigned long iPoint) const override { return ThermalCond(iPoint); }

  /*!
   * \brief Get the vib-el. thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity_ve(unsigned long iPoint) const override { return ThermalCond_ve(iPoint); }

  /*!
   * \brief Set the temperature at the wall
   */
  inline void SetWallTemperature(unsigned long iPoint, su2double temperature_wall) override {
    Primitive(iPoint,T_INDEX) = temperature_wall;
  }

  /*!
   * \brief Get the value of the vorticity.
   * \return Value of the vorticity.
   */
  inline su2double *GetVorticity(unsigned long iPoint) override { return Vorticity[iPoint]; }

};
