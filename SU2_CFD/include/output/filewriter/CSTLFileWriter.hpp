/*!
 * \file CSTLFileWriter.hpp
 * \brief Headers fo the STL file writer class.
 * \author T. Kattmann, T. Albring
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

#include "CFileWriter.hpp"

/*!
 * \class CSTLFileWriter
 * \brief Class for writing STL output files.
 * \author T. Kattmann, T. Albring
 * \version 8.0.0 "Harrier"
 */
class CSTLFileWriter final : public CFileWriter{
private:

  /*--- Variables for reprocessing element connectivity. ---*/

  vector<unsigned long> halo_nodes; /*!< \brief Set of global renumbered node numbers which the current rank refers to via an element, but does not own. Happens on process borders.
                                                Only one process has the 'border' element so that it is not printed twice. */
  unsigned long num_halo_nodes; /*!< \brief Total number of halo nodes this process should receive, i.e. number it needs to acces but doesn't know the data of. */
  /*--- Fields for the first MPI::Alltoallv: Afterwards each process knows the global renumbered node ID's it has to send and to which proc to send. ---*/
  vector<unsigned long> sorted_halo_nodes; /*!< \brief Holds same data as `halo_nodes` but as a vector instead of a set. */
  vector<int> num_nodes_to_receive; /*!< \brief Number of points to receive from each process. */
  vector<int> nodes_to_receive_displacements; /*!< \brief Displacements vector for MPI::AlltoAllv command. */
  vector<unsigned long> nodes_to_send; /*!< \brief Vector of global renumbered node numbers which the current rank owns, and another rank needs for connectivity. */
  vector<int> num_nodes_to_send; /*!< \brief Number of points to to send to each process. */
  vector<int> nodes_to_send_displacements; /*!< \brief Displacements vector for MPI::AlltoAllv command. */
  /*--- Fields for the second MPI::Alltoallv: Afterwards, `halo_var_data` holds the correct coordinate data of halo points of this process. ---*/
  vector<passivedouble> data_to_send; /*!< \brief Holds the halo-point coordinates for other ranks which are owned by this rank. These are communicated. */
  vector<int> num_values_to_send; /*!< \brief Number of Coord values to send to each process. */
  vector<int> values_to_send_displacements; /*!< \brief Displacements vector for MPI::AlltoAllv command. */
  vector<passivedouble> halo_var_data; /*!< \brief Holds the halo-point coordinates for this rank received from other ranks. */
  vector<int> num_values_to_receive; /*!< \brief Number of Coord values to receive from each process. */
  vector<int> values_to_receive_displacements; /*!< \brief Displacements vector for MPI::AlltoAllv command. */

  /*--- Variables for gathering the triangle data in one array. ---*/
  unsigned long *buffRecvTriaCount = nullptr; /*!< \brief Array with number of triangles which each processor has. (Note: Quads are split into two Tris)  */
  unsigned long max_nLocalTriaAll; /*!< \brief Largest Tri count of all processors. */
  su2double *buffSendCoords = nullptr; /*!< \brief Array holding Coordinate data of one processor. 3 consecutive doubles make a point and 3 consecutive points make a Tri. */
  su2double *buffRecvCoords = nullptr; /*!< \brief Array holding Coordinate data of all processors. 3 consecutive doubles make a point and 3 consecutive points make a Tri. */

public:

  const static string fileExt; /*!< \brief File extension ".stl". */

private:

  /*!
   * \brief Recompute Tri/Quad element connectivity between processor borders.
   */
  void ReprocessElementConnectivity();

  /*!
   * \brief Create an array containing all coordinate data for the surface triangles.
   */
  void GatherCoordData();

  /*!
   * \brief Write element coordinate data into a send-buffer in 'triangle order' (i.e. 3 points with 3 coords consecutively).
   * \param[in] elemType - Either TRIANGLE or QUADRILATERAL
   * \param[in] nLocalElements - number of elements of elemType
   * \param[in] startIndex - index to start writing into buffSendCoords
   */
  void StoreCoordData(enum GEO_TYPE elemType,
                      unsigned long nLocalElements,
                      unsigned long startIndex);

  /*!
   * \brief Get the halo-node value of a global renumbered Point for a specific variable.
   * \param[in] global_node_number - global renumbered Point ID
   * \param[in] iVar - Variable number
   * \return Value of the halo-node variable
   */
  passivedouble GetHaloNodeValue(unsigned long global_node_number, unsigned short iVar);

public:

  /*!
   * \brief Construct a file writer using field name and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CSTLFileWriter(CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CSTLFileWriter() override;

  /*!
   * \brief Write sorted data to file in STL file format
   * \param[in] val_filename - The name of the file
   */
  void WriteData(string val_filename) override ;

};

