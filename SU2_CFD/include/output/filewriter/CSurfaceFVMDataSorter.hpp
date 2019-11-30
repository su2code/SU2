/*!
 * \file CSurfaceFVMDataSorter.hpp
 * \brief Headers for the surface FVM data sorter class.
 * \author T. Albring
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

#include "CParallelDataSorter.hpp"
#include "CFVMDataSorter.hpp"

class CSurfaceFVMDataSorter final: public CParallelDataSorter{

  CFVMDataSorter* volume_sorter;                    //!< Pointer to the volume sorter instance
  //! Structure to map the local sorted point ID to the global point ID
  map<unsigned long,unsigned long> Renumber2Global;
public:

  /*!
   * \brief Constructor
   * \param[in] config - Pointer to the current config structure
   * \param[in] geometry - Pointer to the current geometry
   * \param[in] nFields - Number of output fields
   * \param[in] volume_sorter - Pointer to the corresponding volume sorter instance
   */
  CSurfaceFVMDataSorter(CConfig *config, CGeometry* geometry, unsigned short nFields, CFVMDataSorter* volume_sorter);

  /*!
   * \brief Destructor
   */
  ~CSurfaceFVMDataSorter() override;

  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData() override;

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] surf - boolean controlling whether surface <TRUE> or volume connectivity <FALSE> should be sorted.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) override;

  /*!
   * \brief Get the global index of a point.
   * \input iPoint - the point ID.
   * \return Global index of a specific point.
   */
  unsigned long GetGlobalIndex(unsigned long iPoint)  override{
    return Renumber2Global[iPoint];
  }


private:

  /*!
   * \brief Sort the connectivity for a single surface element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

};
