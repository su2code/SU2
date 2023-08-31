/*!
 * \file CIncEulerVariable.hpp
 * \brief Class for defining the variables of the incompressible Euler solver.
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

#include <limits>
#include "CFlowVariable.hpp"

/*!
 * \class CIncEulerVariable
 * \brief Class for defining the variables of the incompressible Euler solver.
 * \note Primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv)
 * \note Gradients of primitives (P, vx, vy, vz, T, rho, beta)
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncEulerVariable : public CFlowVariable {
public:
  static constexpr size_t MAXNVAR = 12;

  template <class IndexType>
  struct CIndices {
    const IndexType nDim;
    CIndices(IndexType ndim, IndexType) : nDim(ndim) {}
    inline IndexType NDim() const { return nDim; }
    inline IndexType NSpecies() const { return 0; }
    inline IndexType Pressure() const { return 0; }
    inline IndexType Velocity() const { return 1; }
    inline IndexType Temperature() const { return nDim+1; }
    inline IndexType Density() const { return nDim+2; }
    inline IndexType Beta() const { return nDim+3; }
    inline IndexType SoundSpeed() const { return Beta(); }
    inline IndexType LaminarViscosity() const { return nDim+4; }
    inline IndexType EddyViscosity() const { return nDim+5; }
    inline IndexType ThermalConductivity() const { return nDim+6; }
    inline IndexType CpTotal() const { return nDim+7; }
    inline IndexType CvTotal() const { return nDim+8; }

    /*--- For compatible interface with NEMO. ---*/
    inline IndexType SpeciesDensities() const { return std::numeric_limits<IndexType>::max(); }
    inline IndexType Temperature_ve() const { return std::numeric_limits<IndexType>::max(); }
    inline IndexType Enthalpy() const { return std::numeric_limits<IndexType>::max(); }
  };

 protected:
  const CIndices<unsigned long> indices;

  VectorType Streamwise_Periodic_RecoveredPressure,    /*!< \brief Recovered/Physical pressure [Pa] for streamwise periodic flow. */
             Streamwise_Periodic_RecoveredTemperature; /*!< \brief Recovered/Physical temperature [K] for streamwise periodic flow. */
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_pressure - value of the pressure.
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] temperature - Value of the temperature (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncEulerVariable(su2double pressure, const su2double *velocity, su2double temperature,
                    unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config);

  /*!
   * \brief Set the value of the pressure.
   * \param[in] iPoint - Point index.
   */
  inline void SetPressure(unsigned long iPoint) final { Primitive(iPoint, indices.Pressure()) = Solution(iPoint,0); }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline bool SetDensity(unsigned long iPoint, su2double val_density) final {
    Primitive(iPoint, indices.Density()) = val_density;
    return val_density <= 0.0;
  }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline void SetVelocity(unsigned long iPoint) final {
    Velocity2(iPoint) = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      Primitive(iPoint, iDim+indices.Velocity()) = Solution(iPoint,iDim+1);
      Velocity2(iPoint) += pow(Primitive(iPoint, iDim+indices.Velocity()), 2);
    }
  }

  /*!
   * \brief Set the value of the temperature for incompressible flows with energy equation.
   * \param[in] iPoint - Point index.
   */
  inline bool SetTemperature(unsigned long iPoint, su2double val_temperature) final {
    Primitive(iPoint, indices.Temperature()) = val_temperature;
    return val_temperature <= 0.0;
  }

  /*!
   * \brief Set the value of the beta coeffient for incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline void SetBetaInc2(unsigned long iPoint, su2double betainc2) final {
    Primitive(iPoint, indices.Beta()) = betainc2;
  }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(unsigned long iPoint) const final { return Primitive(iPoint, indices.Pressure()); }

  /*!
   * \brief Get the value of beta squared for the incompressible flow
   * \return Value of beta squared.
   */
  inline su2double GetBetaInc2(unsigned long iPoint) const final { return Primitive(iPoint, indices.Beta()); }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(unsigned long iPoint) const final { return Primitive(iPoint, indices.Density()); }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(unsigned long iPoint) const final { return Primitive(iPoint, indices.Temperature()); }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final {
    return Primitive(iPoint, iDim+indices.Velocity());
  }

  /*!
   * \brief Get the velocity gradient.
   * \return Value of the velocity gradient.
   */
  inline CMatrixView<const su2double> GetVelocityGradient(unsigned long iPoint) const final {
    return Gradient_Primitive(iPoint, indices.Velocity());
  }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const final {
    su2double ProjVel = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      ProjVel += Primitive(iPoint, iDim+indices.Velocity())*val_vector[iDim];
    return ProjVel;
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Old(iPoint,iDim+1) = val_velocity[iDim];
  }

  /*!
   * \brief Set the momentum part of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVel_ResTruncError_Zero(unsigned long iPoint) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Res_TruncError(iPoint,iDim+1) = 0.0;
  }

  /*!
   * \brief Set all the primitive variables for incompressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) final;

  /*!
   * \brief Set the specific heat Cp.
   */
  inline void SetSpecificHeatCp(unsigned long iPoint, su2double val_Cp) final {
    Primitive(iPoint, indices.CpTotal()) = val_Cp;
  }

  /*!
   * \brief Set the specific heat Cv.
   */
  inline void SetSpecificHeatCv(unsigned long iPoint, su2double val_Cv) final {
    Primitive(iPoint, indices.CvTotal()) = val_Cv;
  }

  /*!
   * \brief Get the specific heat at constant P of the flow.
   * \return Value of the specific heat at constant P of the flow.
   */
  inline su2double GetSpecificHeatCp(unsigned long iPoint) const final { return Primitive(iPoint, indices.CpTotal()); }

  /*!
   * \brief Get the specific heat at constant V of the flow.
   * \return Value of the specific heat at constant V of the flow.
   */
  inline su2double GetSpecificHeatCv(unsigned long iPoint) const final { return Primitive(iPoint, indices.CvTotal()); }

  /*!
   * \brief Set the recovered pressure for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \param[in] val_pressure - pressure value.
   */
  inline void SetStreamwise_Periodic_RecoveredPressure(unsigned long iPoint, su2double val_pressure) final {
    Streamwise_Periodic_RecoveredPressure(iPoint) = val_pressure;
  }

  /*!
   * \brief Get the recovered pressure for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \return Recovered/Physical pressure for streamwise periodic flow.
   */
  inline su2double GetStreamwise_Periodic_RecoveredPressure(unsigned long iPoint) const final {
    return Streamwise_Periodic_RecoveredPressure(iPoint);
  }

  /*!
   * \brief Set the recovered temperature for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \param[in] val_temperature - temperature value.
   */
  inline void SetStreamwise_Periodic_RecoveredTemperature(unsigned long iPoint, su2double val_temperature) final {
    Streamwise_Periodic_RecoveredTemperature(iPoint) = val_temperature;
  }

  /*!
   * \brief Get the recovered temperature for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \return Recovered/Physical temperature for streamwise periodic flow.
   */
  inline su2double GetStreamwise_Periodic_RecoveredTemperature(unsigned long iPoint) const final {
    return Streamwise_Periodic_RecoveredTemperature(iPoint);
  }

  /*!
   * \brief Specify a vector to set the velocity components of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionVector(unsigned long iPoint, const su2double *val_vector) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Solution(iPoint, iDim+1) = val_vector[iDim];
  }

};
