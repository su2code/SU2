/*!
 * \file CIncNSVariable.hpp
 * \brief Class for defining the variables of the incompressible
          Navier-Stokes solver.
 * \author F. Palacios, T. Economon
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

#include "CIncEulerVariable.hpp"

/*!
 * \class CIncNSVariable
 * \brief Class for defining the variables of the incompressible Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncNSVariable : public CIncEulerVariable {
private:
  su2double Vorticity[3];    /*!< \brief Vorticity of the fluid. */
  su2double StrainMag;       /*!< \brief Magnitude of rate of strain tensor. */

  su2double DES_LengthScale;
public:

  /*!
   * \brief Constructor of the class.
   */
  CIncNSVariable(void);

  /*!
   * \overload
   * \param[in] val_pressure - value of the pressure.
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_temperature - Value of the temperature (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncNSVariable(su2double val_pressure, su2double *val_velocity, su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncNSVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CIncNSVariable(void);

  /*!
   * \brief Set the laminar viscosity.
   */
  inline void SetLaminarViscosity(su2double laminarViscosity) {Primitive[nDim+4] = laminarViscosity;}

  /*!
   * \brief Set the vorticity value.
   */
  bool SetVorticity(void);

  /*!
   * \brief Set the rate of strain magnitude.
   */
  bool SetStrainMag(void);

  /*!
   * \overload
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline void SetEddyViscosity(su2double eddy_visc) {Primitive[nDim+5] = eddy_visc; }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(void) {return Primitive[nDim+4]; }

  /*!
   * \brief Get the eddy viscosity of the flow.
   * \return The eddy viscosity of the flow.
   */
  inline su2double GetEddyViscosity(void) {return Primitive[nDim+5]; }

  /*!
   * \brief Set the thermal conductivity.
   */
  inline void SetThermalConductivity(su2double thermalConductivity) {Primitive[nDim+6] = thermalConductivity;}

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(void) {return Primitive[nDim+6]; }

  /*!
   * \brief Get the value of the vorticity.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the vorticity.
   */
  inline su2double *GetVorticity(void) {return Vorticity; }

  /*!
   * \brief Get the value of the magnitude of rate of strain.
   * \return Value of the rate of strain magnitude.
   */
  inline su2double GetStrainMag(void) {return StrainMag; }

  /*!
   * \brief Set all the primitive variables for incompressible flows
   */
  bool SetPrimVar(su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel);
  using CVariable::SetPrimVar;

  /*!
   * \brief Set the DES Length Scale.
   */
  inline void SetDES_LengthScale(su2double val_des_lengthscale) {DES_LengthScale = val_des_lengthscale; }

  /*!
   * \brief Get the DES length scale
   * \return Value of the DES length Scale.
   */
  inline su2double GetDES_LengthScale(void) {return DES_LengthScale; }

};
