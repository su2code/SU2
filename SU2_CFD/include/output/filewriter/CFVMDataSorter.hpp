/*!
 * \file CFVMDataSorter.hpp
 * \brief Headers fo the FVM data sorter class.
 * \author T. Albring
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

#include "CParallelDataSorter.hpp"
#include <vector>

class CFVMDataSorter final: public CParallelDataSorter{
  
private:
  
  int* Local_Halo; //!< Array containing the flag whether a point is a halo node

public:
  
  /*!
   * \brief Constructor
   * \param[in] config - Pointer to the current config structure
   * \param[in] geometry - Pointer to the current geometry
   * \param[in] nFields - Number of output fields
   */
  CFVMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields);
  
  /*!
   * \brief Destructor
   */
  ~CFVMDataSorter() override;
  
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
  unsigned long GetGlobalIndex(unsigned long iPoint) override { 
    return linearPartitioner->GetFirstIndexOnRank(rank) + iPoint;
  }
  
  /*!
   * \brief Get the boolean whether a point is a halo node
   * \param[in] iPoint - ID of the point
   * \return <TRUE> if the point is a halo node.
   */
  bool GetHalo(unsigned long iPoint) {return Local_Halo[iPoint];}
  
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
