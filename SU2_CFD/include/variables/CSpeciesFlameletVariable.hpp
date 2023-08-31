/*!
 * \file CSpeciesFlameletVariable.hpp
 * \brief Base class for defining the variables of the flamelet transport model.
 * \author D. Mayer, T. Economon, N. Beishuizen
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

#include "CSpeciesVariable.hpp"

/*!
 * \class CSpeciesFlameletVariable
 * \brief Base class for defining the variables of the flamelet model.
 */
class CSpeciesFlameletVariable final : public CSpeciesVariable {
 protected:
  MatrixType source_scalar; /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  MatrixType lookup_scalar; /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  su2vector<unsigned short> table_misses; /*!< \brief Vector of lookup table misses. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] species_inf - species variable values (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSpeciesFlameletVariable(const su2double* species_inf, unsigned long npoint, unsigned long ndim, unsigned long nvar,
                           const CConfig* config);

  /*!
   * \brief Set the value of the transported scalar source term.
   * \param[in] iPoint - the location where the value has to be set.
   * \param[in] val_lookup_scalar - the value of the scalar to set.
   * \param[in] val_ivar - Eqn. index to the transport equation.
   */
  inline void SetLookupScalar(unsigned long iPoint, su2double val_lookup_scalar, unsigned short val_ivar) override {
    lookup_scalar(iPoint, val_ivar) = val_lookup_scalar;
  }

  /*!

   * \brief Set a source term for the specie transport equation.
   * \param[in] iPoint - Node index.
   * \param[in] val_ivar - Species index.
   * \param[in] val_source - Source term value.
  */
  inline void SetScalarSource(unsigned long iPoint, unsigned short val_ivar, su2double val_source) override {
    source_scalar(iPoint, val_ivar) = val_source;
  }

  /*!
   * \brief Get the value of the transported scalars source term.
   * \return Pointer to the transported scalars source term.
   */
  inline const su2double* GetScalarSources(unsigned long iPoint) const override { return source_scalar[iPoint]; }

  /*!
   * \brief Get the value of the looked up table based on the transported scalar.
   * \return Pointer to the transported scalars source term.
   */
  inline const su2double* GetScalarLookups(unsigned long iPoint) const override { return lookup_scalar[iPoint]; }

  inline void SetTableMisses(unsigned long iPoint, unsigned short misses) override { table_misses[iPoint] = misses; }

  inline unsigned short GetTableMisses(unsigned long iPoint) const override { return table_misses[iPoint]; }
};
