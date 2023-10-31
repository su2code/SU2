/*!
 * \file CNSVariable.hpp
 * \brief Class for defining the variables of the compressible Navier-Stokes solver.
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

#include "CEulerVariable.hpp"

/*!
 * \class CNSVariable
 * \brief Class for defining the variables of the compressible Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios, T. Economon
 */
class CNSVariable final : public CEulerVariable {
private:
  su2double inv_TimeScale;   /*!< \brief Inverse of the reference time scale. */

  VectorType Tau_Wall;        /*!< \brief Magnitude of the wall shear stress from a wall function. */
  VectorType DES_LengthScale; /*!< \brief DES Length Scale. */
  VectorType Roe_Dissipation; /*!< \brief Roe low dissipation coefficient. */
  VectorType Vortex_Tilting;  /*!< \brief Value of the vortex tilting variable for DES length scale computation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] density - Value of the flow density (initialization value).
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] energy - Value of the flow energy (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNSVariable(su2double density, const su2double *velocity, su2double energy,
              unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config);

  /*!
   * \brief Set the laminar viscosity.
   */
  inline void SetLaminarViscosity(unsigned long iPoint, su2double laminarViscosity) override {
    Primitive(iPoint, indices.LaminarViscosity()) = laminarViscosity;
  }

  /*!
   * \brief Set the laminar viscosity.
   */
  inline void SetThermalConductivity(unsigned long iPoint, su2double thermalConductivity) override {
    Primitive(iPoint, indices.ThermalConductivity()) = thermalConductivity;
  }

  /*!
   * \brief Set the specific heat Cp.
   */
  inline void SetSpecificHeatCp(unsigned long iPoint, su2double val_Cp) override {
    Primitive(iPoint, indices.CpTotal()) = val_Cp;
  }

  /*!
   * \overload
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline void SetEddyViscosity(unsigned long iPoint, su2double eddy_visc) override {
    Primitive(iPoint, indices.EddyViscosity()) = eddy_visc;
  }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.LaminarViscosity());
  }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.ThermalConductivity());
  }

  /*!
   * \brief Get the eddy viscosity of the flow.
   * \return The eddy viscosity of the flow.
   */
  inline su2double GetEddyViscosity(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.EddyViscosity());
  }

  /*!
   * \brief Get the specific heat at constant P of the flow.
   * \return Value of the specific heat at constant P  of the flow.
   */
  inline su2double GetSpecificHeatCp(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.CpTotal());
  }

  /*!
   * \brief Set the derivative of temperature with respect to density (at constant internal energy).
   */
  inline void SetdTdrho_e(unsigned long iPoint, su2double dTdrho_e) override { Secondary(iPoint,2) = dTdrho_e;}

  /*!
   * \brief Set the derivative of temperature with respect to internal energy (at constant density).
   */
  inline void SetdTde_rho(unsigned long iPoint, su2double dTde_rho) override { Secondary(iPoint,3) = dTde_rho;}

  /*!
   * \brief Set the derivative of laminar viscosity with respect to density (at constant temperature).
   */
  inline void Setdmudrho_T(unsigned long iPoint, su2double dmudrho_T) override { Secondary(iPoint,4) = dmudrho_T;}

  /*!
   * \brief Set the derivative of laminar viscosity with respect to temperature (at constant density).
   */
  inline void SetdmudT_rho(unsigned long iPoint, su2double dmudT_rho) override { Secondary(iPoint,5) = dmudT_rho;}

  /*!
   * \brief Set the derivative of thermal conductivity with respect to density (at constant temperature).
   */
  inline void Setdktdrho_T(unsigned long iPoint, su2double dktdrho_T) override { Secondary(iPoint,6) = dktdrho_T;}

  /*!
   * \brief Set the derivative of thermal conductivity with respect to temperature (at constant density).
   */
  inline void SetdktdT_rho(unsigned long iPoint, su2double dktdT_rho) override { Secondary(iPoint,7) = dktdT_rho;}

  /*!
   * \brief Set all the primitive variables for compressible flows
   */
  bool SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) override;
  using CVariable::SetPrimVar;

  /*!
   * \brief Set all the secondary variables (partial derivatives) for compressible flows
   */
  void SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) override;

  /*!
   * \brief Set the value of the wall shear stress computed by a wall function.
   */
  inline void SetTau_Wall(unsigned long iPoint, su2double tau_wall) override { Tau_Wall(iPoint) = tau_wall; }

  /*!
   * \brief Get the value of the wall shear stress computed by a wall function.
   * \return Value of the wall shear stress computed by a wall function.
   */
  inline su2double GetTau_Wall(unsigned long iPoint) const override { return Tau_Wall(iPoint); }
  inline const VectorType& GetTau_Wall() const  { return Tau_Wall; }

  /*!
   * \brief Get the DES length scale
   * \return Value of the DES length Scale.
   */
  inline su2double GetDES_LengthScale(unsigned long iPoint) const override { return DES_LengthScale(iPoint); }

  /*!
   * \brief Set the DES Length Scale.
   */
  inline void SetDES_LengthScale(unsigned long iPoint, su2double val_des_lengthscale) override {
    DES_LengthScale(iPoint) = val_des_lengthscale;
  }

  /*!
   * \brief Set the new solution for Roe Dissipation.
   * \param[in] val_delta - A scalar measure of the grid size
   * \param[in] val_const_DES - The DES constant (C_DES)
   */
  void SetRoe_Dissipation_NTS(unsigned long iPoint, su2double val_delta, su2double val_const_DES) override;

  /*!
   * \brief Set the new solution for Roe Dissipation.
   */
  void SetRoe_Dissipation_FD(unsigned long iPoint, su2double wall_distance) override;

  /*!
   * \brief Get the Roe Dissipation Coefficient.
   * \return Value of the Roe Dissipation.
   */
  inline su2double GetRoe_Dissipation(unsigned long iPoint) const override { return Roe_Dissipation(iPoint); }
  inline const VectorType& GetRoe_Dissipation() const { return Roe_Dissipation; }

  /*!
   * \brief Set the Roe Dissipation Coefficient.
   * \param[in] val_dissipation - Value of the Roe dissipation factor.
   */
  inline void SetRoe_Dissipation(unsigned long iPoint, su2double val_dissipation) override {
    Roe_Dissipation(iPoint) = val_dissipation;
  }

};
