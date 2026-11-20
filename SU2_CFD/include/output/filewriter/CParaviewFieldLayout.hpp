/*!
 * \file CParaviewFieldLayout.hpp
 * \brief Selection of the scalar and vector fields written by the Paraview file writers.
 * \author A. Rausa
 * \version 8.6.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include <algorithm>
#include <string>
#include <vector>

/*!
 * \struct CParaviewField
 * \brief A field as written by the Paraview writers: a scalar or a vector of 3 components.
 */
struct CParaviewField {
  std::size_t first;         /*!< \brief Index of the first output field of this field. */
  unsigned short components; /*!< \brief 1 for a scalar, 3 for a vector (z is 0 in 2D). */
  std::string name;          /*!< \brief Name without quotes, and without "_x" for a vector. */
};

/*!
 * \brief Group the output fields into the scalars and vectors written by the Paraview writers.
 * \note The coordinates (the first nDim fields) are skipped. A field "A_x" followed by "A_y" (and "A_z" in 3D)
 *       is a vector "A". Any other field, also a vector component selected on its own, is a scalar.
 * \param[in] fieldNames - Names of the output fields.
 * \param[in] nDim - Number of dimensions.
 * \returns The fields to write, in the order of the output fields.
 */
inline std::vector<CParaviewField> GetParaviewFields(const std::vector<std::string>& fieldNames,
                                                     unsigned short nDim) {
  auto names = fieldNames;
  for (auto& name : names) name.erase(std::remove(name.begin(), name.end(), '"'), name.end());

  std::vector<CParaviewField> fields;
  for (std::size_t i = nDim; i < names.size();) {
    const auto& name = names[i];
    const bool startsVector = name.size() > 2 && name.compare(name.size() - 2, 2, "_x") == 0;
    const auto base = startsVector ? name.substr(0, name.size() - 2) : name;
    const bool isVector = (nDim == 2 || nDim == 3) && startsVector && i + nDim <= names.size() &&
                          names[i + 1] == base + "_y" && (nDim == 2 || names[i + 2] == base + "_z");
    fields.push_back({i, static_cast<unsigned short>(isVector ? 3 : 1), isVector ? base : name});
    i += isVector ? nDim : 1;
  }
  return fields;
}
