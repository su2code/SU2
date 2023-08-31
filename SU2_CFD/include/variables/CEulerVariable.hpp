/*!
 * \file CEulerVariable.hpp
 * \brief Class for defining the variables of the compressible Euler solver.
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
 * \class CEulerVariable
 * \brief Class for defining the variables of the compressible Euler solver.
 * \note Primitive variables (T, vx, vy, vz, P, rho, h, c)
 * \note Gradients and limiters (T, vx, vy, vz, P, rho)
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon
 */
class CEulerVariable : public CFlowVariable {
 public:
  static constexpr size_t MAXNVAR = 12;

  template <class IndexType>
  struct CIndices {
    const IndexType nDim;
    CIndices(IndexType ndim, IndexType) : nDim(ndim) {}
    inline IndexType NDim() const { return nDim; }
    inline IndexType NSpecies() const { return 0; }
    inline IndexType Temperature() const { return 0; }
    inline IndexType Velocity() const { return 1; }
    inline IndexType Pressure() const { return nDim+1; }
    inline IndexType Density() const { return nDim+2; }
    inline IndexType Enthalpy() const { return nDim+3; }
    inline IndexType SoundSpeed() const { return nDim+4; }
    inline IndexType LaminarViscosity() const { return nDim+5; }
    inline IndexType EddyViscosity() const { return nDim+6; }
    inline IndexType ThermalConductivity() const { return nDim+7; }
    inline IndexType CpTotal() const { return nDim+8; }

    /*--- For compatible interface with NEMO. ---*/
    inline IndexType SpeciesDensities() const { return std::numeric_limits<IndexType>::max(); }
    inline IndexType Temperature_ve() const { return std::numeric_limits<IndexType>::max(); }
  };

 protected:
  const CIndices<unsigned long> indices;

  /*!< \brief Secondary variables (dPdrho_e, dPde_rho, dTdrho_e, dTde_rho, dmudrho_T, dmudT_rho, dktdrho_T, dktdT_rho)
   *          in compressible (Euler: 2, NS: 8) flows. */
  MatrixType Secondary;

  MatrixType WindGust;      /*! < \brief Wind gust value */

  bool DataDrivenFluid = false; /*!< \brief Usage of data-driven fluid model. DatasetExtrapolation and FluidEntropy will not be sized if disabled. */
  su2vector<unsigned short> DatasetExtrapolation; /*!< \brief Stores instances of dataset bounds violation when using data-driven fluid models. */
  su2vector<unsigned long> NIterNewtonsolver;    /*!< \brief Stores number of Newton solver iterations when using data-driven fluid models. */
  VectorType FluidEntropy;          /*!< \brief Stores the fluid entropy value as computed by the data-driven fluid model. */

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
  CEulerVariable(su2double density, const su2double *velocity, su2double energy,
                 unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config);

  /*!
   * \brief A virtual member.
   */
  inline void SetdPdrho_e(unsigned long iPoint, su2double dPdrho_e) final { Secondary(iPoint,0) = dPdrho_e;}

  /*!
   * \brief A virtual member.
   */
  inline void SetdPde_rho(unsigned long iPoint, su2double dPde_rho) final { Secondary(iPoint,1) = dPde_rho;}

  /*!
   * \brief Set the value of the pressure.
   */
  inline bool SetPressure(unsigned long iPoint, su2double pressure) final {
    Primitive(iPoint,nDim+1) = pressure;
    return pressure <= 0.0;
  }

  /*!
   * \brief Set the value of the speed of the sound.
   * \param[in] soundspeed2 - Value of soundspeed^2.
   */
  bool SetSoundSpeed(unsigned long iPoint, su2double soundspeed2) final {
    if (soundspeed2 < 0.0) return true;
    else {
      Primitive(iPoint,nDim+4) = sqrt(soundspeed2);
      return false;
    }
  }

  /*!
   * \brief Set the value of the enthalpy.
   */
  inline void SetEnthalpy(unsigned long iPoint) final {
    Primitive(iPoint, indices.Enthalpy()) =
      (Solution(iPoint,nVar-1) + Primitive(iPoint, indices.Pressure())) / Solution(iPoint,0);
  }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) final;

  /*!
   * \brief A virtual member.
   */
  void SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) override;

  /*!
   * \brief Get all the secondary variables.
   */
  inline const MatrixType& GetSecondary() const {return Secondary; }

  /*!
   * \brief Get the secondary variables.
   * \param[in] iVar - Index of the variable.
   * \return Value of the secondary variable for the index <i>iVar</i>.
   */
  inline su2double GetSecondary(unsigned long iPoint, unsigned long iVar) const final { return Secondary(iPoint,iVar); }

  /*!
   * \brief Set the value of the secondary variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the secondary variable for the index <i>iVar</i>.
   */
  inline void SetSecondary(unsigned long iPoint, unsigned long iVar, su2double val_secondary) final {
    Secondary(iPoint,iVar) = val_secondary;
  }

  /*!
   * \brief Set the value of the secondary variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the secondary variable for the index <i>iVar</i>.
   */
  inline void SetSecondary(unsigned long iPoint, const su2double *val_secondary) final {
    for (unsigned long iVar = 0; iVar < nSecondaryVar; iVar++)
      Secondary(iPoint,iVar) = val_secondary[iVar];
  }

  /*!
   * \brief Get the secondary variables of the problem.
   * \return Pointer to the secondary variable vector.
   */
  inline su2double *GetSecondary(unsigned long iPoint) final { return Secondary[iPoint]; }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   */
  inline bool SetDensity(unsigned long iPoint) final {
    Primitive(iPoint, indices.Density()) = Solution(iPoint,0);
    return Primitive(iPoint, indices.Density()) <= 0.0;
  }

  /*!
   * \brief Set the value of the temperature.
   * \param[in] temperature - how agitated the particles are :)
   */
  inline bool SetTemperature(unsigned long iPoint, su2double temperature) final {
    Primitive(iPoint, indices.Temperature()) = temperature;
    return temperature <= 0.0;
  }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(unsigned long iPoint) const final { return Primitive(iPoint, indices.Pressure()); }

  /*!
   * \brief Get the speed of the sound.
   * \return Value of speed of the sound.
   */
  inline su2double GetSoundSpeed(unsigned long iPoint) const final { return Primitive(iPoint, indices.SoundSpeed()); }

  /*!
   * \brief Get the enthalpy of the flow.
   * \return Value of the enthalpy of the flow.
   */
  inline su2double GetEnthalpy(unsigned long iPoint) const final { return Primitive(iPoint, indices.Enthalpy()); }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(unsigned long iPoint) const final { return Solution(iPoint,0); }

  /*!
   * \brief Get the energy of the flow.
   * \return Value of the energy of the flow.
   */
  inline su2double GetEnergy(unsigned long iPoint) const final { return Solution(iPoint,nVar-1)/Solution(iPoint,0); }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(unsigned long iPoint) const final { return Primitive(iPoint,indices.Temperature()); }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final {
    return Primitive(iPoint,iDim+indices.Velocity());
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
      ProjVel += Primitive(iPoint,iDim+indices.Velocity())*val_vector[iDim];
    return ProjVel;
  }

  /*!
   * \brief Set the velocity vector from the solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity(unsigned long iPoint) final {
    Velocity2(iPoint) = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      Primitive(iPoint,iDim+indices.Velocity()) = Solution(iPoint,iDim+1) / Solution(iPoint,0);
      Velocity2(iPoint) += pow(Primitive(iPoint,iDim+indices.Velocity()),2);
    }
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Old(iPoint,iDim+1) = val_velocity[iDim]*Solution(iPoint,0);
  }

  /*!
   * \brief Set the momentum part of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVel_ResTruncError_Zero(unsigned long iPoint) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Res_TruncError(iPoint,iDim+1) = 0.0;
  }

  /*!
   * \brief Specify a vector to set the velocity components of the solution. Multiplied by density for compressible cases.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionVector(unsigned long iPoint, const su2double *val_vector) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Solution(iPoint, iDim+1) = GetDensity(iPoint) * val_vector[iDim];
  }

  /*!
   * \brief Set fluid entropy
   * \param[in] iPoint - Node index
   * \param[in] entropy - fluid entropy value.
   */
  inline void SetEntropy(unsigned long iPoint, su2double entropy) final { FluidEntropy[iPoint] = entropy; };

  /*!
   * \brief Get fluid entropy
   * \param[in] iPoint - Node index
   * \return Entropy - Fluid entropy value
   */
  inline su2double GetEntropy(unsigned long iPoint) const final { return FluidEntropy[iPoint]; }

  /*!
   * \brief Set dataset extrapolation instance
   * \param[in] iPoint - Node index
   * \param[in] extrapolation - Extrapolation instance (0 = within dataset, 1 = outside dataset)
   */
  inline void SetDataExtrapolation(unsigned long iPoint, unsigned short extrapolation) final {
    DatasetExtrapolation[iPoint] = extrapolation;
  };

  /*!
   * \brief Get dataset extrapolation instance
   * \param[in] iPoint - Node index
   * \return extrapolation - Extrapolation instance (0 = within dataset, 1 = outside dataset)
   */
  inline unsigned short GetDataExtrapolation(unsigned long iPoint) const final { return DatasetExtrapolation[iPoint]; }

  /*!
   * \brief Set the number of iterations required by a Newton solver used by the fluid model.
   * \param[in] iPoint - Node index
   * \param[in] nIter - Number of iterations evaluated by the Newton solver
   */
  inline void SetNewtonSolverIterations(unsigned long iPoint, unsigned long nIter) final { NIterNewtonsolver[iPoint] = nIter; }

  /*!
   * \brief Get the number of iterations required by a Newton solver used by the fluid model.
   * \param[in] iPoint - Node index
   * \return Number of iterations evaluated by the Newton solver
   */
  inline unsigned long GetNewtonSolverIterations(unsigned long iPoint) const final { return NIterNewtonsolver[iPoint]; }

};
