/*!
 * \file CSurfaceFVMDataSorter.hpp
 * \brief Headers for the surface FVM data sorter class.
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
#include "CFVMDataSorter.hpp"

class CSurfaceFVMDataSorter final: public CParallelDataSorter{

  const CFVMDataSorter* volumeSorter;               //!< Pointer to the volume sorter instance
  map<unsigned long,unsigned long> Renumber2Global; //! Structure to map the local sorted point ID to the global point ID
public:

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] config - Pointer to the current config structure
   * \param[in] geometry - Pointer to the current geometry
   * \param[in] valVolumeSorter - The datasorter containing the volume data
   */
  CSurfaceFVMDataSorter(CConfig *config, CGeometry* geometry, const CFVMDataSorter* valVolumeSorter);

  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData() override;

  /*!
   * \brief Sort the connectivities on the surface into data structures used for output file writing.
   *  All markers in MARKER_PLOTTING will be sorted.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_sort - boolean controlling whether surface <TRUE> or volume connectivity <FALSE> should be sorted.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) override;

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing.
   * Only markers in the markerList argument will be sorted.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] markerList - List of markers to sort.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, const vector<string> &markerList) override;

  /*!
   * \brief Get the global index of a point.
   * \param[in] iPoint - the local renumbered point ID.
   * \return Global index of a specific point.
   */
  unsigned long GetGlobalIndex(unsigned long iPoint) const override{
    if (iPoint > nPoints)
      SU2_MPI::Error(string("Local renumbered iPoint ID ") + to_string(iPoint) +
                     string(" is larger than max number of nodes ") + to_string(nPoints), CURRENT_FUNCTION);

    return Renumber2Global.at(iPoint);
  }

  /*!
   * \brief Get the beginning global renumbered node ID of the linear partition owned by a specific processor.
   * \param[in] rank - the processor rank.
   * \return The beginning global renumbered node ID.
   */
  unsigned long GetNodeBegin(unsigned short rank) const override {
    return nPoint_Recv[rank];
  }

  /*!
   * \brief Get the Processor ID a Point belongs to.
   * \param[in] iPoint - global renumbered ID of the point
   * \return The rank/processor number.
   */
  unsigned short FindProcessor(unsigned long iPoint) const override {

    for (unsigned short iRank = 1; iRank < size; iRank++){
      if (nPoint_Recv[iRank] > static_cast<int>(iPoint)){
        return iRank - 1;
      }
    }
    return size-1;
  }

  /*!
   * \brief Get the cumulated number of points
   * \input rank - the processor rank.
   * \return The cumulated number of points up to certain processor rank.
   */
  unsigned long GetnPointCumulative(unsigned short rank) const override {
    return GetNodeBegin(rank);
  }

private:

  /*!
   * \brief Sort the connectivity for a single surface element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   * \param[in] markerList - List of markers to sort
   */
  void SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type,
                               const vector<string> &markerList);

};
