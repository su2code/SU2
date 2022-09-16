/*!
 * \file CCGNSMeshReaderFVM.hpp
 * \brief Header file for the class CCGNSMeshReaderFVM.
 *        The implementations are in the <i>CCGNSMeshReaderFVM.cpp</i> file.
 * \author T. Economon
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CCGNSMeshReaderFVM
 * \brief Reads a CGNS zone into linear partitions for the finite volume solver (FVM).
 * \author: T. Economon
 */
class CCGNSMeshReaderFVM final: public CCGNSMeshReaderBase {
  
private:

#ifdef HAVE_CGNS
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

  /*!
   * \brief Routine to launch non-blocking sends and recvs amongst all processors.
   * \param[in] bufSend - Buffer of data to be sent.
   * \param[in] nElemSend - Array containing the number of elements to send to other processors in cumulative storage format.
   * \param[in] sendReq - Array of MPI send requests.
   * \param[in] bufRecv - Buffer of data to be received.
   * \param[in] nElemSend - Array containing the number of elements to receive from other processors in cumulative storage format.
   * \param[in] sendReq - Array of MPI recv requests.
   * \param[in] countPerElem - Pieces of data per element communicated.
   */
  void InitiateCommsAll(void *bufSend,
                        const int *nElemSend,
                        SU2_MPI::Request *sendReq,
                        void *bufRecv,
                        const int *nElemRecv,
                        SU2_MPI::Request *recvReq,
                        unsigned short countPerElem,
                        unsigned short commType);

  /*!
   * \brief Routine to complete the set of non-blocking communications launched with InitiateComms() with MPI_Waitany().
   * \param[in] nSends - Number of sends to be completed.
   * \param[in] sendReq - Array of MPI send requests.
   * \param[in] nRecvs - Number of receives to be completed.
   * \param[in] sendReq - Array of MPI recv requests.
   */
  void CompleteCommsAll(int nSends,
                        SU2_MPI::Request *sendReq,
                        int nRecvs,
                        SU2_MPI::Request *recvReq);


public:

  /*!
   * \brief Constructor of the CCGNSMeshReaderFVM class.
   */
  CCGNSMeshReaderFVM(CConfig        *val_config,
                     unsigned short val_iZone,
                     unsigned short val_nZone);

};
