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
class CTNE2NSVariable : public CTNE2EulerVariable {
private:
  su2double Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
  su2double Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
  su2double Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
  su2double Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
  su2double *DiffusionCoeff;   /*!< \brief Diffusion coefficient of the mixture. */
  su2double **Dij;             /*!< \brief Binary diffusion coefficients. */
  su2double LaminarViscosity;  /*!< \brief Viscosity of the fluid. */
  su2double ThermalCond;       /*!< \brief T-R thermal conductivity of the gas mixture. */
  su2double ThermalCond_ve;    /*!< \brief V-E thermal conductivity of the gas mixture. */
  su2double Vorticity[3];	   /*!< \brief Vorticity of the fluid. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CTNE2NSVariable(void);

  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nVarPrim - Number of primitive variables.
   * \param[in] val_nVarPrimGrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2NSVariable(unsigned short val_nDim, unsigned short val_nVar,
                    unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                    CConfig *config);

  /*!
   * \overload
   * \param[in] val_density - Value of the flow density (initialization value).
   * \param[in] val_massfrac - Value of the flow mass fraction (initialization value).
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_temperature - Value of the flow temperature (initialization value).
   * \param[in] val_temperature_ve - Value of the flow temperature_ve (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nPrimVar - Number of primitive variables.
   * \param[in] val_nPrimVargrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2NSVariable(su2double val_density, su2double *val_massfrac, su2double *val_velocity,
                  su2double val_temperature, su2double val_temperature_ve, unsigned short val_nDim,
                  unsigned short val_nVar, unsigned short val_nPrimVar,
                  unsigned short val_nPrimVarGrad, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nPrimVar - Number of primitive variables.
   * \param[in] val_nPrimgVarGrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2NSVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nVar,
                  unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                  CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTNE2NSVariable(void);

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
  inline su2double* GetDiffusionCoeff(void) { return DiffusionCoeff; }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(void) { return LaminarViscosity; }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(void) { return ThermalCond; }

  /*!
   * \brief Get the vib-el. thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity_ve(void) { return ThermalCond_ve; }

  /*!
   * \brief Set the temperature at the wall
   */
  inline void SetWallTemperature(su2double Temperature_Wall) { Primitive[T_INDEX] = Temperature_Wall; }

  /*!
   * \brief Get the value of the vorticity.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the vorticity.
   */
  inline su2double *GetVorticity(unsigned short val_dim) { return Vorticity; }

  /*!
   * \brief Set all the primitive variables for compressible flows
   */
  bool SetPrimVar_Compressible(CConfig *config);

};
