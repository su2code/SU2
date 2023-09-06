/*!
 * \file CFVMDataSorter.hpp
 * \brief Headers fo the FVM data sorter class.
 * \author T. Albring
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

#include "CParallelDataSorter.hpp"
#include <vector>

class CFVMDataSorter final: public CParallelDataSorter{

private:

  vector<int> Local_Halo; //!< Array containing the flag whether a point is a halo node

public:
  /*!
   * \brief Constructor
   * \param[in] config - Pointer to the current config structure
   * \param[in] geometry - Pointer to the current geometry
   * \param[in] valFieldNames - Vector containing the field names
   */
  CFVMDataSorter(CConfig *config, CGeometry *geometry, const vector<string> &valFieldNames);

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) override;

  /*!
   * \brief Get the global index of a point.
   * \input iPoint - the point ID.
   * \return Global index of a specific point.
   */
  unsigned long GetGlobalIndex(unsigned long iPoint) const override {
    return linearPartitioner.GetFirstIndexOnRank(rank) + iPoint;
  }

  /*!
   * \brief Get the boolean whether a point is a halo node
   * \param[in] iPoint - ID of the point
   * \return <TRUE> if the point is a halo node.
   */
  bool GetHalo(unsigned long iPoint) const {return Local_Halo[iPoint];}

private:

  /*!
   * \brief Initialize the halo point flags
   * \param[in] geometry - Pointer to the current geometry
   * \param[in] config - Pointer to the current config structure
   */
  void SetHaloPoints(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Sort the connectivity for a single volume element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  void SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type, bool val_sort);

};
