/*!
 * \file CParallelDataSorter.hpp
 * \brief Headers fo the data sorter class.
 * \author T. Albring, T. Economon
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

#include "../../../../Common/include/mpi_structure.hpp"
#include "../../../../Common/include/option_structure.hpp"
#include "../../../../Common/include/toolboxes/CLinearPartitioner.hpp"

class CGeometry;
class CConfig;

class CParallelDataSorter{
protected:

  /*!
   * \brief The MPI rank
   */
  int rank;

  /*!
   * \brief The MPI size, aka the number of processors.
   */
  int size;

  unsigned long nGlobal_Poin_Par;  //!< Global number of points without halos before sorting
  unsigned long nGlobal_Elem_Par;  //!< Global number of elems without halos before sorting
  unsigned long nParallel_Poin;    //!< Local number of points after sorting on this proc
  unsigned long nParallel_Line,    //!< Local number of line elements after sorting on this proc
  nParallel_Tria,                  //!< Local number of triangle elements after sorting on this proc
  nParallel_Quad,                  //!< Local number of quad elements after sorting on this proc
  nParallel_Tetr,                  //!< Local number of tetrahedral elements after sorting on this proc
  nParallel_Hexa,                  //!< Local number of hexhedral elements after sorting on this proc
  nParallel_Pris,                  //!< Local number of prism elements after sorting on this proc
  nParallel_Pyra;                  //!< Local number of pyramid elements after sorting on this proc
  int *Conn_Line_Par;              //!< Local connectivity of line elements after sorting on this proc
  int *Conn_Tria_Par;              //!< Local connectivity of triangle elements after sorting on this proc
  int *Conn_Quad_Par;              //!< Local connectivity of quad elements after sorting on this proc
  int *Conn_Tetr_Par;              //!< Local connectivity of tetrahedral elements after sorting on this proc
  int *Conn_Hexa_Par;              //!< Local connectivity of hexahedral elements after sorting on this proc
  int *Conn_Pris_Par;              //!< Local connectivity of prism elements after sorting on this proc
  int *Conn_Pyra_Par;              //!< Local connectivity of pyramid elements after sorting on this proc

  unsigned long nGlobalPoint_Sort; //!< Global number of points without halos after sorting
  unsigned long nLocalPoint_Sort;  //!< Local number of points without halos after sorting on this proc


  CLinearPartitioner* linearPartitioner;  //!< Linear partitioner based on the global number of points.

  unsigned short GlobalField_Counter;  //!< Number of output fields

  bool connectivity_sorted;            //!< Boolean to store information on whether the connectivity is sorted

  int *nPoint_Send;                    //!< Number of points this processor has to send to other processors
  int *nPoint_Recv;                    //!< Number of points this processor receives from other processors
  unsigned long *Index;                //!< Index each point has in the send buffer
  su2double *connSend;                 //!< Send buffer holding the data that will be send to other processors
  passivedouble *passiveDoubleBuffer;  //!< Buffer holding the sorted, partitioned data as passivedouble types
  su2double     *doubleBuffer;         //!< Buffer holding the sorted, partitioned data as su2double types
  /// Pointer used to allocate the memory used for ::passiveDoubleBuffer and ::doubleBuffer.
  char *dataBuffer;
  unsigned long *idSend;               //!< Send buffer holding global indices that will be send to other processors
  int nSends,                          //!< Number of sends
  nRecvs;                              //!< Number of receives

  /*!
   * \brief Prepare the send buffers by filling them with the global indices.
   * After calling this function, the data buffer for sending can be filled with the
   * ::SetUnsorted_Data() routine.
   * \param[in] globalID - Vector containing the global indices of the points
   */
  void PrepareSendBuffers(std::vector<unsigned long>& globalID);

public:

  /*!
   * \brief Constructor
   * \param[in] config - Pointer to the current config structure
   * \param[in] nFields - Number of output fields
   */
  CParallelDataSorter(CConfig *config, unsigned short nFields);

  /*!
   * \brief Destructor
   */
  virtual ~CParallelDataSorter();

  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors.
   */
  virtual void SortOutputData();

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  virtual void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort = true){}

  /*!
   * \brief Get the number of points the local rank owns.
   * \return local number of points.
   */
  unsigned long GetnPoints(){return nParallel_Poin;}

  /*!
   * \brief Get the number of points to sort.
   * \return local number of points.
   */
  unsigned long GetnLocalPointSort(){return nLocalPoint_Sort;}

  /*!
   * \brief Get the global number of points (accumulated from all ranks)
   * \return Global number of points.
   */
  unsigned long GetnPointsGlobal(){return nGlobal_Poin_Par;}

  /*!
   * \brief Get the global of elements (accumulated from all ranks and element types)
   * \return Global number elements.
   */
  unsigned long GetnElem(){return nGlobal_Elem_Par;}

  /*!
   * \brief Get the local number of elements of a specific type that the current rank owns
   * \input type - The type of element, ref GEO_TYPE
   * \return Local number of elements of a specific type.
   */
  unsigned long GetnElem(GEO_TYPE type);

  /*!
   * \brief Get the connectivity of specific element.
   * \input type - The type of element, ref GEO_TYPE
   * \input iElem - The element ID
   * \input iNode - The node ID
   * \return the connected node.
   */
  unsigned long GetElem_Connectivity(GEO_TYPE type, unsigned long iElem, unsigned long iNode);

  /*!
   * \brief Beginning node ID of the linear partition owned by a specific processor.
   * \input rank - the processor rank.
   * \return The beginning node ID.
   */
  unsigned long GetNodeBegin(unsigned short rank){return linearPartitioner->GetFirstIndexOnRank(rank);}

  /*!
   * \brief Ending node ID of the linear partition owned by a specific processor.
   * \input rank - the processor rank.
   * \return The ending node ID.
   */
  unsigned long GetNodeEnd(unsigned short rank){return linearPartitioner->GetLastIndexOnRank(rank);}

  /*!
   * \brief Get the value of the linear partitioned data.
   * \input iField - the output field ID.
   * \input iPoint - the point ID.
   * \return the value of the data field at a point.
   */
  passivedouble GetData(unsigned short iField, unsigned long iPoint) {return passiveDoubleBuffer[iPoint*GlobalField_Counter + iField];}

  /*!
   * \brief Get the pointer to the sorted linear partitioned data.
   * \return Pointer to the sorted data.
   */
  const passivedouble *GetData() {return passiveDoubleBuffer;}

  /*!
   * \brief Get the global index of a point.
   * \input iPoint - the point ID.
   * \return Global index of a specific point.
   */
  virtual unsigned long GetGlobalIndex(unsigned long iPoint){return 0;}

  /*!
   * \brief Get the cumulated number of points
   * \input rank - the processor rank.
   * \return The cumulated number of points up to certain processor rank.
   */
  unsigned long GetnPointCumulative(unsigned short rank){return linearPartitioner->GetCumulativeSizeBeforeRank(rank);}

  /*!
   * \brief Get the linear number of points
   * \input rank - the processor rank.
   * \return The linear number of points up to certain processor rank.
   */
  unsigned long GetnPointLinear(unsigned short rank){return linearPartitioner->GetSizeOnRank(rank);}

  /*!
   * \brief Check whether the current connectivity is sorted (i.e. if SortConnectivity has been called)
   * \return <TRUE> if the connectivity is sorted.
   */
  bool GetConnectivitySorted(){return connectivity_sorted;}

  /*!
   * \brief Set the value of a specific field at a point.
   * ::PrepareSendBuffers must be called before using this function.
   *
   * \param[in] iPoint - ID of the point
   * \param[in] iField - Index of the field
   * \param[in] data - Value of the field
   */
  void SetUnsorted_Data(unsigned long iPoint, unsigned short iField, su2double data){
    connSend[Index[iPoint] + iField] = data;
  }

  su2double GetUnsorted_Data(unsigned long iPoint, unsigned short iField){
    return connSend[Index[iPoint] + iField];
  }
};
