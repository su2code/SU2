/*!
 * \file CMarkerProfileReaderFVM.hpp
 * \brief Header file for the class CMarkerProfileReaderFVM.
 *        The implementations are in the <i>CMarkerProfileReaderFVM.cpp</i> file.
 * \author T. Economon
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

#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/geometry/CGeometry.hpp"

using namespace std;

/*!
 * \class CMarkerProfileReaderFVM
 * \brief Class for the marker profile reader of the finite volume solver (FVM).
 * \author: T. Economon
 */
class CMarkerProfileReaderFVM {

protected:

  int rank;  /*!< \brief MPI Rank. */
  int size;  /*!< \brief MPI Size. */

  CConfig *config;  /*!< \brief Local pointer to the config parameter object. */

  CGeometry *geometry;  /*!< \brief Local pointer to the geometry object. */

  unsigned short dimension;  /*!< \brief Dimension of the problem (2 or 3). */
  unsigned short markerType;  /*!< \brief Type of marker where the profiles are being applied. */

  unsigned short numberOfVars;  /*!< \brief Number of variables added to the number of coordinates to write each line in the template profile file. */

  unsigned long numberOfProfiles;  /*!< \brief Auxiliary structure for holding the number of markers in a profile file. */

  string filename;  /*!< \brief File name of the marker profile file. */

  vector<string> columnNames; /*!< \brief string containing all the names of the columns, one for each marker */
  vector<string> columnValues; /*!< \brief initial values for the profile, in string format */

  vector<string> profileTags;  /*!< \brief Auxiliary structure for holding the string names of the markers in a profile file. */

  vector<unsigned long> numberOfRowsInProfile;  /*!< \brief Auxiliary structure for holding the number of rows for a particular marker in a profile file. */
  vector<unsigned long> numberOfColumnsInProfile;  /*!< \brief Auxiliary structure for holding the number of columns for a particular marker in a profile file. */
  vector<string> totalColumnNames; /*!< \brief Names of the columns for the profile, one for each inlet marker. */
  vector<string> totalColumnValues; /*!< \brief Initial values for the profile, constructed from MARKER_INLET. */

  vector<vector<passivedouble> > profileData;  /*!< \brief Auxiliary structure for holding the data values from a profile file. */
  vector<vector<vector<su2double> > > profileCoords;  /*!< \brief Data structure for holding the merged inlet boundary coordinates from all ranks. */

private:

  /*!
   * \brief Read a native SU2 marker profile file in ASCII format.
   */
  void ReadMarkerProfile();

  /*!
   * \brief Merge the node coordinates of all profile-type boundaries from all processors.
   */
  void MergeProfileMarkers();

  /*!
   * \brief Write a template profile file if the requested file is not found.
   */
  void WriteMarkerProfileTemplate();

public:

  /*!
   * \brief Constructor of the CMarkerProfileReaderFVM class.
   * \param[in] val_geometry    - Pointer to the current geometry
   * \param[in] val_config      - Pointer to the current config structure
   * \param[in] val_filename    - Name of the profile file to be read
   * \param[in] val_kind_marker - Type of marker where profile will be applied
   * \param[in] val_number_vars - Number of columns of profile data to be written to template file (excluding coordinates)
   */
  CMarkerProfileReaderFVM(CGeometry      *val_geometry,
                          CConfig        *val_config,
                          string         val_filename,
                          unsigned short val_kind_marker,
                          unsigned short val_number_vars,
                          vector<string> val_columnNames,
                          vector<string> val_columnValues);

  /*!
   * \brief Destructor of the CMeshReaderFVM class.
   */
  ~CMarkerProfileReaderFVM(void);

  /*!
   * \brief Get the number of profiles found within the input file.
   * \returns Number of profiles found within the input file.
   */
  inline unsigned long GetNumberOfProfiles() const {
    return numberOfProfiles;
  }

  /*!
   * \brief Get the string tag for the marker where the profile is applied.
   * \param[in] val_iProfile - current profile index.
   * \returns String tag for the marker where the profile is applied.
   */
  inline const string &GetTagForProfile(int val_iProfile) const {
    return profileTags[val_iProfile];
  }

  /*!
   * \brief Get the number of rows of data in a profile.
   * \param[in] val_iProfile - current profile index.
   * \returns Number of rows of data in a profile.
   */
  inline unsigned long GetNumberOfRowsInProfile(int val_iProfile) {
    return numberOfRowsInProfile[val_iProfile];
  }

  /*!
   * \brief Get the number of columns of data in a profile.
   * \param[in] val_iProfile - current profile index.
   * \returns Number of columns of data in a profile.
   */
  inline unsigned long GetNumberOfColumnsInProfile(int val_iProfile) {
    return numberOfColumnsInProfile[val_iProfile];
  }

  /*!
   * \brief Get the 1D vector of data for a profile from the input file.
   * \param[in] val_iProfile - current profile index.
   * \returns 1D vector of data for a profile from the input file.
   */
  inline const vector<passivedouble> &GetDataForProfile(int val_iProfile) const {
    return profileData[val_iProfile];
  }

   /*!
   * \brief Get the data for the specific column if interpolation being done.
   * \param[in] val_iProfile - current profile index.
   * \param[in] iCol - the column whose data is required
   * \returns the specific column data.
   */
  inline vector<su2double> GetColumnForProfile(int val_iProfile, unsigned short iCol) const {
    auto nRow = numberOfRowsInProfile[val_iProfile];
    auto nCol = numberOfColumnsInProfile[val_iProfile];
    vector<su2double> ColumnData(nRow);
    for (unsigned long iRow = 0; iRow < nRow; iRow++)
      ColumnData[iRow]=profileData[val_iProfile][iRow*nCol+iCol];
    return ColumnData;
  }
};
