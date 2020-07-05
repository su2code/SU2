/*!
 * \file CCGNSMeshReaderFEM.hpp
 * \brief Header file for the class CCGNSMeshReaderFEM.
 *        The implementations are in the <i>CCGNSMeshReaderFEM.cpp</i> file.
 * \author T. Economon, E. van der Weide
 * \version 7.0.5 "Blackbird"
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

#include "CCGNSMeshReaderBase.hpp"

/*!
 * \class CCGNSMeshReaderFEM
 * \brief Reads a CGNS zone into linear partitions for the finite element solver (FEM).
 * \author: T. Economon, E. van der Weide
 */
class CCGNSMeshReaderFEM: public CCGNSMeshReaderBase {
  
private:
  
#ifdef HAVE_CGNS
  
  /*!
   * \brief Reads the grid points from a CGNS zone into linear partitions across all ranks.
   */
  void ReadCGNSPointCoordinates();
  
  /*!
   * \brief Reads the interior volume elements from one section of a CGNS zone into linear partitions across all ranks.
   * \param[in] val_section - CGNS section index.
   */
  void ReadCGNSVolumeSection(int val_section);
  
  /*!
   * \brief Reads the surface (boundary) elements from the CGNS zone. Only the master rank currently reads and stores the connectivity, which is linearly partitioned later.
   * \param[in] val_section - CGNS section index.
   */
  void ReadCGNSSurfaceSection(int val_section);
  
  /*!
   * \brief Reformats the CGNS volume connectivity from file into the standard base class data structures.
   */
  void ReformatCGNSVolumeConnectivity();
  
  /*!
   * \brief Reformats the CGNS volume connectivity from file into the standard base class data structures.
   */
  void ReformatCGNSSurfaceConnectivity();
  
#endif
  
public:
  
  /*!
   * \brief Constructor of the CCGNSMeshReaderFEM class.
   */
  CCGNSMeshReaderFEM(CConfig        *val_config,
                     unsigned short val_iZone,
                     unsigned short val_nZone);
  
  /*!
   * \brief Destructor of the CCGNSMeshReaderFEM class.
   */
  ~CCGNSMeshReaderFEM(void);
  
};
