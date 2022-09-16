/*!
 * \file CCGNSMeshReaderBase.hpp
 * \brief Header file for the class CCGNSMeshReaderBase.
 *        The implementations are in the <i>CCGNSMeshReaderBase.cpp</i> file.
 * \author T. Economon
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#ifdef HAVE_CGNS
#include "cgnslib.h"
#endif

#include "CMeshReaderFVM.hpp"

/*!
 * \class CCGNSMeshReaderBase
 * \brief Base class for the reading of a CGNS zone.
 * \author: T. Economon
 */
class CCGNSMeshReaderBase: public CMeshReaderFVM {
  
protected:
 
#ifdef HAVE_CGNS
  int cgnsFileID; /*!< \brief CGNS file identifier. */
  int cgnsBase;   /*!< \brief CGNS database index. */
  int cgnsZone;   /*!< \brief CGNS zone index. */
  
  int nZones;     /*!< \brief Total number of zones in the CGNS file. */
  int nSections;  /*!< \brief Total number of sections in the CGNS file. */
  
  vector<bool> isInterior;             /*!< \brief Vector of booleans to store whether each section in the CGNS file is an interior or boundary section. */
  vector<unsigned long> nElems;        /*!< \brief Vector containing the local number of elements found within each CGNS section. */
  vector<unsigned long> elemOffset;    /*!< \brief Global ID offset for each interior section (i.e., the total number of global elements that came before it). */
  vector<vector<cgsize_t> > connElems; /*!< \brief Vector containing the local element connectivity found within each CGNS section. First index is the section, second contains the connectivity in format [globalID VTK n1 n2 n3 n4 n5 n6 n7 n8] for each element. */
  vector<vector<char> > sectionNames;  /*!< \brief Vector for storing the names of each boundary section (marker). */
  
  /*!
   * \brief Open the CGNS file and checks for errors.
   * \param[in] val_filename - string name of the CGNS file to be read.
   */
  void OpenCGNSFile(string val_filename);
  
  /*!
   * \brief Reads all CGNS database metadata and checks for errors.
   */
  void ReadCGNSDatabaseMetadata();

  /*!
   * \brief Reads the grid points from a CGNS zone into linear partitions across all ranks.
   */
  void ReadCGNSPointCoordinates();

  /*!
   * \brief Reads all CGNS zone metadata and checks for errors.
   */
  void ReadCGNSZoneMetadata();

  /*!
   * \brief Reads the metadata for each CGNS section in a zone and collect information, including the size and whether it is an interior or boundary section.
   */
  void ReadCGNSSectionMetadata();
  
  /*!
   * \brief Get the VTK type and string name for a CGNS element type.
   * \param[in] val_elem_type - CGNS element type.
   * \param[out] val_vtk_type - VTK type identifier index.
   * \returns String containing the name of the element type.
   */
  string GetCGNSElementType(ElementType_t val_elem_type,
                            int           &val_vtk_type);
#endif
  
public:
  
  /*!
   * \brief Constructor of the CCGNSMeshReaderBase class.
   */
  CCGNSMeshReaderBase(CConfig        *val_config,
                      unsigned short val_iZone,
                      unsigned short val_nZone);
  
  /*!
   * \brief Destructor of the CCGNSMeshReaderBase class.
   */
  virtual ~CCGNSMeshReaderBase(void);
  
};
