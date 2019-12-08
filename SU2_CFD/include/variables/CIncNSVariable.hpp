/*!
 * \file CIncNSVariable.hpp
 * \brief Class for defining the variables of the incompressible
          Navier-Stokes solver.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
class CIncNSVariable final : public CIncEulerVariable {
private:
  MatrixType Vorticity;    /*!< \brief Vorticity of the fluid. */
  VectorType StrainMag;    /*!< \brief Magnitude of rate of strain tensor. */

  VectorType DES_LengthScale;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] pressure - value of the pressure.
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] temperature - Value of the temperature (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncNSVariable(su2double pressure, const su2double *velocity, su2double temperature,
                 unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CIncNSVariable() = default;

  /*!
   * \brief Set the laminar viscosity.
   */
  inline void SetLaminarViscosity(unsigned long iPoint, su2double laminarViscosity) override {
    Primitive(iPoint,nDim+4) = laminarViscosity;
  }

  /*!
   * \brief Set the vorticity value.
   */
  bool SetVorticity_StrainMag() override;

  /*!
   * \overload
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline void SetEddyViscosity(unsigned long iPoint, su2double eddy_visc) override {
    Primitive(iPoint,nDim+5) = eddy_visc;
  }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(unsigned long iPoint) const override { return Primitive(iPoint,nDim+4); }

  /*!
   * \brief Get the eddy viscosity of the flow.
   * \return The eddy viscosity of the flow.
   */
  inline su2double GetEddyViscosity(unsigned long iPoint) const override { return Primitive(iPoint,nDim+5); }

  /*!
   * \brief Set the thermal conductivity.
   */
  inline void SetThermalConductivity(unsigned long iPoint, su2double thermalConductivity) override {
    Primitive(iPoint,nDim+6) = thermalConductivity;
  }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(unsigned long iPoint) const override { return Primitive(iPoint,nDim+6); }

  /*!
   * \brief Get the value of the vorticity.
   * \return Value of the vorticity.
   */
  inline su2double *GetVorticity(unsigned long iPoint) override { return Vorticity[iPoint]; }

  /*!
   * \brief Get the value of the magnitude of rate of strain.
   * \return Value of the rate of strain magnitude.
   */
  inline su2double GetStrainMag(unsigned long iPoint) const override { return StrainMag(iPoint); }

  /*!
   * \brief Set all the primitive variables for incompressible flows
   */
  bool SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) override;
  using CVariable::SetPrimVar;

  /*!
   * \brief Set the DES Length Scale.
   */
  inline void SetDES_LengthScale(unsigned long iPoint, su2double val_des_lengthscale) override {
    DES_LengthScale(iPoint) = val_des_lengthscale;
  }

  /*!
   * \brief Get the DES length scale
   * \return Value of the DES length Scale.
   */
  inline su2double GetDES_LengthScale(unsigned long iPoint) const override { return DES_LengthScale(iPoint); }

};
