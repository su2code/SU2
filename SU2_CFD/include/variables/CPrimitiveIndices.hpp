/*!
 * \file CPrimitiveIndices.hpp
 * \brief Abstract representation of flow primitive variable indices that tries to be efficient.
 * \version 7.5.0 "Blackbird"
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

#include <cstdint>
#include <cstring>
#include "CEulerVariable.hpp"
#include "CIncEulerVariable.hpp"
#include "CNEMOEulerVariable.hpp"

/*!
 * \class CPrimitiveIndices
 * \brief Abstract wrapper for the CIndices classes of CEulerVariable, CIncEulerVariable, etc.
 * This should only be used when the concrete type of solver variable is not known.
 * \note We are using some tricks to try to make this more efficient than virtual functions.
 */
template <class IndexType>
struct CPrimitiveIndices {
  CPrimitiveIndices(bool incompressible, bool nemo, IndexType nDim, IndexType nSpecies) {
    /*--- Instantiate the right type of CIndices, without heap allocations. ---*/
    if (incompressible) {
      type_ = 0;
      Construct<CIncEulerVariable::template CIndices<IndexType>>(nDim, nSpecies);
    } else if (nemo) {
      type_ = 1;
      Construct<CNEMOEulerVariable::template CIndices<IndexType>>(nDim, nSpecies);
    } else {
      type_ = 2;
      Construct<CEulerVariable::template CIndices<IndexType>>(nDim, nSpecies);
    }
  }

  /*--- To use the CIndices object we stored on the buffer, we cast to the
   * concrete type using the same logic used during construction. ---*/
#define GET_INDEX_IMPL(NAME)                                                                            \
  inline IndexType NAME() const {                                                                       \
    switch (type_) {                                                                                    \
     case 0:                                                                                            \
      return reinterpret_cast<const CIncEulerVariable::template CIndices<IndexType>*>(data_)->NAME();   \
     case 1:                                                                                            \
      return reinterpret_cast<const CNEMOEulerVariable::template CIndices<IndexType>*>(data_)->NAME();  \
     default:                                                                                           \
      return reinterpret_cast<const CEulerVariable::template CIndices<IndexType>*>(data_)->NAME();      \
    }                                                                                                   \
  }
  GET_INDEX_IMPL(Temperature)
  GET_INDEX_IMPL(Velocity)
  GET_INDEX_IMPL(Pressure)
  GET_INDEX_IMPL(Density)
  GET_INDEX_IMPL(Enthalpy)
  GET_INDEX_IMPL(SoundSpeed)
  GET_INDEX_IMPL(LaminarViscosity)
  GET_INDEX_IMPL(EddyViscosity)
  GET_INDEX_IMPL(ThermalConductivity)
  GET_INDEX_IMPL(CpTotal)
  GET_INDEX_IMPL(Temperature_ve)

 private:
  template <class ConcreteIndices>
  void Construct(IndexType nDim, IndexType nSpecies) {
    /*--- Build the indices object in the static buffer owned by this class. ---*/
    static_assert(sizeof(ConcreteIndices) <= 2 * sizeof(IndexType), "");
    new(data_) ConcreteIndices(nDim, nSpecies);
  }

  alignas(sizeof(IndexType)) char data_[2 * sizeof(IndexType)]{};  /*!< \brief Space for the largest CIndices class. */
  uint8_t type_;
};