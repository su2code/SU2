/*!
 * \file CSU2ASCIIMeshReaderFVM.hpp
 * \brief Header file for the class CSU2ASCIIMeshReaderFVM.
 *        The implementations are in the <i>CSU2ASCIIMeshReaderFVM.cpp</i> file.
 * \author T. Economon
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

#include "CMeshReaderFVM.hpp"

/*!
 * \class CSU2ASCIIMeshReaderFVM
 * \brief Reads a native SU2 ASCII grid into linear partitions for the finite volume solver (FVM).
 * \author: T. Economon
 */
class CSU2ASCIIMeshReaderFVM: public CMeshReaderFVM {
  
private:
  
  unsigned short myZone; /*!< \brief Current SU2 zone index. */
  unsigned short nZones; /*!< \brief Total number of zones in the SU2 file. */
  
  string meshFilename; /*!< \brief Name of the SU2 ASCII mesh file being read. */
  ifstream mesh_file;  /*!< \brief File object for the SU2 ASCII mesh file. */
  
  bool actuator_disk; /*!< \brief Boolean for whether we have an actuator disk to split. */
  
  unsigned long ActDiskNewPoints; /*!< \brief Total number of new grid points to add due to actuator disk splitting. */
  
  su2double Xloc; /*!< \brief X-coordinate of the CG of the actuator disk surface. */
  su2double Yloc; /*!< \brief X-coordinate of the CG of the actuator disk surface. */
  su2double Zloc; /*!< \brief X-coordinate of the CG of the actuator disk surface. */
  
  vector<bool> ActDisk_Bool; /*!< \brief Flag to identify the grid points on the actuator disk. */
  
  vector<unsigned long> ActDiskPoint_Back; /*!< \brief Vector containing the global index for the new grid points added to the back of the actuator disk. */
  vector<unsigned long> VolumePoint_Inv; /*!< \brief Vector containing the inverse mapping from the global index to the added point index for the actuator disk. */
  
  vector<su2double> CoordXActDisk; /*!< \brief X-coordinates of the new grid points added by splitting the actuator disk (size = ActDiskNewPoints). */
  vector<su2double> CoordYActDisk; /*!< \brief Y-coordinates of the new grid points added by splitting the actuator disk (size = ActDiskNewPoints). */
  vector<su2double> CoordZActDisk; /*!< \brief Z-coordinates of the new grid points added by splitting the actuator disk (size = ActDiskNewPoints). */
  
  vector<su2double> CoordXVolumePoint; /*!< \brief X-coordinates of the volume elements touching the actuator disk. */
  vector<su2double> CoordYVolumePoint; /*!< \brief Y-coordinates of the volume elements touching the actuator disk. */
  vector<su2double> CoordZVolumePoint; /*!< \brief Z-coordinates of the volume elements touching the actuator disk. */
  
  /*!
   * \brief Reads all SU2 ASCII mesh metadata and checks for errors.
   */
  void ReadMetadata();
  
  /*!
   * \brief Splits a single surface actuator disk boundary into two separate markers (repeated points).
   */
  void SplitActuatorDiskSurface();
  
  /*!
   * \brief Reads the grid points from an SU2 zone into linear partitions across all ranks.
   */
  void ReadPointCoordinates();
  
  /*!
   * \brief Reads the interior volume elements from one section of an SU2 zone into linear partitions across all ranks.
   */
  void ReadVolumeElementConnectivity();
  
  /*!
   * \brief Reads the surface (boundary) elements from the SU2 zone.
   */
  void ReadSurfaceElementConnectivity();
  
  /*!
   * \brief Helper function to find the current zone in an SU2 ASCII mesh object.
   */
  void FastForwardToMyZone();
  
public:
  
  /*!
   * \brief Constructor of the CSU2ASCIIMeshReaderFVM class.
   */
  CSU2ASCIIMeshReaderFVM(CConfig        *val_config,
                         unsigned short val_iZone,
                         unsigned short val_nZone);
  
  /*!
   * \brief Destructor of the CSU2ASCIIMeshReaderFVM class.
   */
  ~CSU2ASCIIMeshReaderFVM(void);
  
};
