/*!
 * \file CCGNSMeshReaderFVM.hpp
 * \brief Header file for the class CCGNSMeshReaderFVM.
 *        The implementations are in the <i>CCGNSMeshReaderFVM.cpp</i> file.
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

#ifdef HAVE_CGNS
#include "cgnslib.h"
#endif

#include "CMeshReaderFVM.hpp"

/*!
 * \class CCGNSMeshReaderFVM
 * \brief Reads a CGNS zone into linear partitions for the finite volume solver (FVM).
 * \author: T. Economon
 */
class CCGNSMeshReaderFVM : public CMeshReaderFVM {
 private:
#ifdef HAVE_CGNS
  int cgnsFileID;         /*!< \brief CGNS file identifier. */
  const int cgnsBase = 1; /*!< \brief CGNS database index (the CGNS reader currently assumes a single database). */
  const int cgnsZone = 1; /*!< \brief CGNS zone index (and 1 zone in that database). */

  int nSections; /*!< \brief Total number of sections in the CGNS file. */

  vector<bool> isInterior; /*!< \brief Vector of booleans to store whether each section in the CGNS file is an interior
                              or boundary section. */
  vector<unsigned long>
      nElems; /*!< \brief Vector containing the local number of elements found within each CGNS section. */
  vector<unsigned long> elemOffset;    /*!< \brief Global ID offset for each interior section (i.e., the total number of
                                          global elements that came before it). */
  vector<vector<cgsize_t> > connElems; /*!< \brief Vector containing the local element connectivity found within each
                                          CGNS section. First index is the section, second contains the connectivity in
                                          format [globalID VTK n1 n2 n3 n4 n5 n6 n7 n8] for each element. */
  vector<vector<char> > sectionNames;  /*!< \brief Vector for storing the names of each boundary section (marker). */

  /*!
   * \brief Open the CGNS file and checks for errors.
   * \param[in] val_filename - string name of the CGNS file to be read.
   */
  void OpenCGNSFile(const string& val_filename);

  /*!
   * \brief Reads all CGNS database metadata and checks for errors.
   */
  void ReadCGNSDatabaseMetadata();

  /*!
   * \brief Reads all CGNS zone metadata and checks for errors.
   */
  void ReadCGNSZoneMetadata();

  /*!
   * \brief Reads the grid points from a CGNS zone into linear partitions across all ranks.
   */
  void ReadCGNSPointCoordinates();

  /*!
   * \brief Reads the metadata for each CGNS section in a zone and collect information, including the size and whether
   * it is an interior or boundary section.
   */
  void ReadCGNSSectionMetadata();

  /*!
   * \brief Reads the interior volume elements from one section of a CGNS zone into linear partitions across all ranks.
   * \param[in] val_section - CGNS section index.
   */
  void ReadCGNSVolumeSection(int val_section);

  /*!
   * \brief Reads the surface (boundary) elements from the CGNS zone. Only the master rank currently reads and stores
   * the connectivity, which is linearly partitioned later. \param[in] val_section - CGNS section index.
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

  /*!
   * \brief Get the VTK type and string name for a CGNS element type.
   * \param[in] val_elem_type - CGNS element type.
   * \param[out] val_vtk_type - VTK type identifier index.
   * \returns String containing the name of the element type.
   */
  string GetCGNSElementType(ElementType_t val_elem_type, int& val_vtk_type);
#endif

  /*!
   * \brief Routine to launch non-blocking sends and recvs amongst all processors.
   * \param[in] bufSend - Buffer of data to be sent.
   * \param[in] nElemSend - Array containing the number of elements to send to other processors in cumulative storage
   * format. \param[in] sendReq - Array of MPI send requests. \param[in] bufRecv - Buffer of data to be received.
   * \param[in] nElemSend - Array containing the number of elements to receive from other processors in cumulative
   * storage format. \param[in] sendReq - Array of MPI recv requests. \param[in] countPerElem - Pieces of data per
   * element communicated.
   */
  void InitiateCommsAll(void* bufSend, const int* nElemSend, SU2_MPI::Request* sendReq, void* bufRecv,
                        const int* nElemRecv, SU2_MPI::Request* recvReq, unsigned short countPerElem,
                        unsigned short commType);

  /*!
   * \brief Routine to complete the set of non-blocking communications launched with InitiateComms() with MPI_Waitany().
   * \param[in] nSends - Number of sends to be completed.
   * \param[in] sendReq - Array of MPI send requests.
   * \param[in] nRecvs - Number of receives to be completed.
   * \param[in] sendReq - Array of MPI recv requests.
   */
  void CompleteCommsAll(int nSends, SU2_MPI::Request* sendReq, int nRecvs, SU2_MPI::Request* recvReq);

 public:
  /*!
   * \brief Constructor of the CCGNSMeshReaderFVM class.
   */
  CCGNSMeshReaderFVM(CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Destructor of the CCGNSMeshReaderFVM class.
   */
  ~CCGNSMeshReaderFVM(void);
};
