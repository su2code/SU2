/*!
 * \file CPrimitiveIndices.hpp
 * \brief Abstract representation of flow primitive variable indices that tries to be efficient.
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

#include <cstdint>
#include <cstring>
#include <map>
#include <string>

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
  GET_INDEX_IMPL(NDim)
  GET_INDEX_IMPL(NSpecies)
  GET_INDEX_IMPL(SpeciesDensities)
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
#undef GET_INDEX_IMPL

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

/*!
 * \brief Maps primitive variable names to their indices based on a CPrimitiveIndices object.
 * \note The variables names we try to map here must be kept in sync with CPrimitiveIndices.
 * All index classes should use numeric_limits::max for variables that are not available.
 */
template <class IndexType>
std::map<std::string, IndexType> PrimitiveNameToIndexMap(const CPrimitiveIndices<IndexType>& idx) {
  std::map<std::string, IndexType> nameToIndex;
  const auto notAvailable = std::numeric_limits<IndexType>::max();
#define INSERT_VAR_INDEX_PAIR(STRING, FUNCTION) \
  if (idx.FUNCTION() != notAvailable) nameToIndex[STRING] = idx.FUNCTION();
  INSERT_VAR_INDEX_PAIR("TEMPERATURE", Temperature)
  INSERT_VAR_INDEX_PAIR("TEMPERATURE_VE", Temperature_ve)
  INSERT_VAR_INDEX_PAIR("PRESSURE", Pressure)
  INSERT_VAR_INDEX_PAIR("DENSITY", Density)
  INSERT_VAR_INDEX_PAIR("ENTHALPY", Enthalpy)
  INSERT_VAR_INDEX_PAIR("SOUND_SPEED", SoundSpeed)
  INSERT_VAR_INDEX_PAIR("LAMINAR_VISCOSITY", LaminarViscosity)
  INSERT_VAR_INDEX_PAIR("EDDY_VISCOSITY", EddyViscosity)
  INSERT_VAR_INDEX_PAIR("THERMAL_CONDUCTIVITY", ThermalConductivity)
  INSERT_VAR_INDEX_PAIR("CP_TOTAL", CpTotal)
#undef INSERT_VAR_INDEX_PAIR
  nameToIndex["VELOCITY_X"] = idx.Velocity();
  nameToIndex["VELOCITY_Y"] = idx.Velocity() + 1;
  if (idx.NDim() == 3) {
    nameToIndex["VELOCITY_Z"] = idx.Velocity() + 2;
  }
  for (IndexType iSpecies = 0; iSpecies < idx.NSpecies(); ++iSpecies) {
    nameToIndex["DENSITY_" + std::to_string(iSpecies)] = idx.SpeciesDensities() + iSpecies;
  }
  return nameToIndex;
}
