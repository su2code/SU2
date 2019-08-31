/*!
 * \file CBoxMeshReaderFVM.hpp
 * \brief Header file for the class CBoxMeshReaderFVM.
 *        The implementations are in the <i>CBoxMeshReaderFVM.cpp</i> file.
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
 * \class CBoxMeshReaderFVM
 * \brief Reads a 3D box grid into linear partitions for the finite volume solver (FVM).
 * \author: T. Economon
 */
class CBoxMeshReaderFVM: public CMeshReaderFVM {
  
private:
  
  unsigned long nNode; /*!< \brief Number of grid nodes in the x-direction. */
  unsigned long mNode; /*!< \brief Number of grid nodes in the y-direction. */
  unsigned long pNode; /*!< \brief Number of grid nodes in the z-direction. */
  
  su2double Lx; /*!< \brief Length of the domain in the x-direction. */
  su2double Ly; /*!< \brief Length of the domain in the y-direction. */
  su2double Lz; /*!< \brief Length of the domain in the z-direction. */
  
  su2double Ox; /*!< \brief Offset of the domain from 0.0 in the x-direction. */
  su2double Oy; /*!< \brief Offset of the domain from 0.0 in the y-direction. */
  su2double Oz; /*!< \brief Offset of the domain from 0.0 in the z-direction. */
  
  unsigned short KindElem;  /*!< \brief VTK identifier of the interior elements. */
  unsigned short KindBound; /*!< \brief VTK identifier of the surface elements. */
  
  /*!
   * \brief Computes and stores the grid points based on an analytic definition of a box grid.
   */
  void ComputeBoxPointCoordinates();
  
  /*!
   * \brief Computes and stores the volume element connectivity based on an analytic definition of a box grid.
   */
  void ComputeBoxVolumeConnectivity();
  
  /*!
   * \brief Computes and stores the surface element connectivity based on an analytic definition of a box grid.
   */
  void ComputeBoxSurfaceConnectivity();
  
public:
  
  /*!
   * \brief Constructor of the CBoxMeshReaderFVM class.
   */
  CBoxMeshReaderFVM(CConfig        *val_config,
                    unsigned short val_iZone,
                    unsigned short val_nZone);
  
  /*!
   * \brief Destructor of the CBoxMeshReaderFVM class.
   */
  ~CBoxMeshReaderFVM(void);
  
};
