/*!
 * \file CParallelDataSorter.hpp
 * \brief Headers fo the data sorter class.
 * \author T. Albring, T. Economon
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

#include "../../../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../../Common/include/option_structure.hpp"
#include "../../../../Common/include/toolboxes/CLinearPartitioner.hpp"
#include <array>
#include <cassert>

class CGeometry;
class CConfig;

class CParallelDataSorter{
protected:

  /*!
   * \brief The MPI rank
   */
  const int rank;

  /*!
   * \brief The MPI size, aka the number of processors.
   */
  const int size;

  unsigned long nGlobalPointBeforeSort; //!< Global number of points without halos before sorting
  unsigned long nLocalPointsBeforeSort;   //!< Local number of points without halos before sorting on this proc

  int *Conn_Line_Par;              //!< Local connectivity of line elements
  int *Conn_Tria_Par;              //!< Local connectivity of triangle elements
  int *Conn_Quad_Par;              //!< Local connectivity of quad elements
  int *Conn_Tetr_Par;              //!< Local connectivity of tetrahedral elements
  int *Conn_Hexa_Par;              //!< Local connectivity of hexahedral elements
  int *Conn_Pris_Par;              //!< Local connectivity of prism elements
  int *Conn_Pyra_Par;              //!< Local connectivity of pyramid elements

  array<unsigned long, N_ELEM_TYPES> nElemPerTypeGlobal;   //!< Global number of elements after sorting on this proc
  array<unsigned long, N_ELEM_TYPES> nElemPerType; //!< Local number of elements after sorting on this proc

  /*!
   * \brief Map that stores the index for each GEO_TYPE type where to find information
   * in the element arrays.
   */
  struct {
    static unsigned short at(unsigned short type) {
      switch(type) {
        case LINE: return 0;
        case TRIANGLE: return 1;
        case QUADRILATERAL: return 2;
        case TETRAHEDRON: return 3;
        case HEXAHEDRON: return 4;
        case PRISM: return 5;
        case PYRAMID: return 6;
        default: assert(false); return 0;
      };
    }
  } TypeMap;

  unsigned long nPointsGlobal;   //!< Global number of points without halos
  unsigned long nElemGlobal;    //!< Global number of elems without halos
  unsigned long nConnGlobal;    //!< Global size of the connectivity array
  unsigned long nPoints;    //!< Local number of points
  unsigned long nElem;     //!< Local number of elements
  unsigned long nConn;     //!< Local size of the connectivity array

  CLinearPartitioner linearPartitioner;  //!< Linear partitioner based on the global number of points.

  unsigned short GlobalField_Counter;  //!< Number of output fields

  bool connectivitySorted;            //!< Boolean to store information on whether the connectivity is sorted

  int *nPoint_Send;                    //!< Number of points this processor has to send to other processors
  int *nPoint_Recv;                    //!< Number of points this processor receives from other processors
  int *nElem_Send;                     //!< Number of elements this processor has to send to other processors
  int *nElem_Cum;                      //!< Cumulative number of elements
  int *nElemConn_Send;                 //!< Number of element connectivity this processor has to send to other processors
  int *nElemConn_Cum;                  //!< Cumulative number of element connectivity entries
  unsigned long *Index;                //!< Index each point has in the send buffer
  passivedouble *connSend;             //!< Send buffer holding the data that will be send to other processors
  passivedouble *dataBuffer;           //!< Buffer holding the sorted, partitioned data as passivedouble types
  unsigned long *idSend;               //!< Send buffer holding global indices that will be send to other processors
  int nSends,                          //!< Number of sends
  nRecvs;                              //!< Number of receives

  vector<string> fieldNames;           //!< Vector with names of the output fields

  unsigned short nDim;                 //!< Spatial dimension of the data

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
   * \param[in] valFieldNames - Vector containing the field names
   */
  CParallelDataSorter(CConfig *config, const vector<string> &valFieldNames);

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
   * \brief Sort the connectivities into data structures (only for surface data sorters).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] markerList - A list of marker names that should be sorted.
   */
  virtual void SortConnectivity(CConfig *config, CGeometry *geometry, const vector<string> &markerList){}

  /*!
   * \brief Get the number of points the local rank owns.
   * \return local number of points.
   */
  unsigned long GetnPoints() const {return nPoints;}

  /*!
   * \brief Get the number of points to sort.
   * \return local number of points.
   */
  unsigned long GetnLocalPointsBeforeSort() const {return nLocalPointsBeforeSort;}

  /*!
   * \brief Get the global number of points (accumulated from all ranks)
   * \return Global number of points.
   */
  unsigned long GetnPointsGlobal() const {return nPointsGlobal;}

  /*!
   * \brief Get the global of elements (accumulated from all ranks and element types)
   * \return Global number elements.
   */
  unsigned long GetnElem() const {return nElem;}

  /*!
   * \brief Get the local number of elements of a specific type that the current rank owns
   * \input type - The type of element, ref GEO_TYPE
   * \return Local number of elements of a specific type.
   */
  unsigned long GetnElem(GEO_TYPE type) const {
    return nElemPerType[TypeMap.at(type)];
  }

  /*!
   * \brief Get the global number of elements of a specific type
   * \input type - The type of element, ref GEO_TYPE
   * \return global number of elements of a specific type.
   */
  unsigned long GetnElemGlobal(GEO_TYPE type) const {
    return nElemPerTypeGlobal[TypeMap.at(type)];
  }

  /*!
   * \brief Get the global number of elements
   * \return global number of elements.
   */
  unsigned long GetnElemGlobal() const {
    return nElemGlobal;
  }

  /*!
   * \brief Get the global number entries of the connectivity array
   * \return global number of entries of the connectivity array
   */
  unsigned long GetnConnGlobal() const {
    return nConnGlobal;
  }

  /*!
   * \brief Get the local number entries of the connectivity array
   * \return local number of entries of the connectivity array
   */
  unsigned long GetnConn() const {
    return nConn;
  }

  /*!
   * \brief Get the cumulated number of elements
   * \input rank - the processor rank.
   * \return The cumulated number of elements
   */
  unsigned long GetnElemCumulative(unsigned short rank) const {
    return nElem_Cum[rank];
  }

  /*!
   * \brief Get the cumulated number of entries of the connectivity array
   * \input rank - the processor rank.
   * \return The cumulated number of entries of the connectivity array
   */
  unsigned long GetnElemConnCumulative(unsigned short rank) const {
    return nElemConn_Cum[rank];
  }

  /*!
   * \brief Get the connectivity of specific element.
   * \input type - The type of element, ref GEO_TYPE
   * \input iElem - The element ID
   * \input iNode - The node ID
   * \return the connected node.
   */
  unsigned long GetElemConnectivity(GEO_TYPE type, unsigned long iElem, unsigned long iNode) const ;

  /*!
   * \brief Beginning node ID of the linear partition owned by a specific processor.
   * \input rank - the processor rank.
   * \return The beginning node ID.
   */
  virtual unsigned long GetNodeBegin(unsigned short rank) const {
    return linearPartitioner.GetFirstIndexOnRank(rank);
  }

  /*!
   * \brief Ending node ID of the linear partition owned by a specific processor.
   * \param rank - the processor rank.
   * \return The ending node ID.
   */
  unsigned long GetNodeEnd(unsigned short rank) const {
    return linearPartitioner.GetLastIndexOnRank(rank);
  }

  /*!
   * \brief Get the value of the linear partitioned data.
   * \input iField - the output field ID.
   * \input iPoint - the point ID.
   * \return the value of the data field at a point.
   */
  passivedouble GetData(unsigned short iField, unsigned long iPoint) const  {return dataBuffer[iPoint*GlobalField_Counter + iField];}

  /*!
   * \brief Get the pointer to the sorted linear partitioned data.
   * \return Pointer to the sorted data.
   */
  const passivedouble *GetData() const {return dataBuffer;}

  /*!
   * \brief Get the global index of a point.
   * \input iPoint - the point ID.
   * \return Global index of a specific point.
   */
  virtual unsigned long GetGlobalIndex(unsigned long iPoint) const { return 0; }

  /*!
   * \brief Get the cumulated number of points
   * \input rank - the processor rank.
   * \return The cumulated number of points up to certain processor rank.
   */
  virtual unsigned long GetnPointCumulative(unsigned short rank) const {return linearPartitioner.GetCumulativeSizeBeforeRank(rank);}

  /*!
   * \brief Get the linear number of points
   * \input rank - the processor rank.
   * \return The linear number of points up to certain processor rank.
   */
  unsigned long GetnPointLinear(unsigned short rank) const {return linearPartitioner.GetSizeOnRank(rank);}

  /*!
   * \brief Check whether the current connectivity is sorted (i.e. if SortConnectivity has been called)
   * \return <TRUE> if the connectivity is sorted.
   */
  bool GetConnectivitySorted() const {return connectivitySorted;}

  /*!
   * \brief Set the value of a specific field at a point.
   * ::PrepareSendBuffers must be called before using this function.
   *
   * \param[in] iPoint - ID of the point
   * \param[in] iField - Index of the field
   * \param[in] data - Value of the field
   */
  void SetUnsortedData(unsigned long iPoint, unsigned short iField, su2double data){
    connSend[Index[iPoint] + iField] = SU2_TYPE::GetValue(data);
  }

  passivedouble GetUnsortedData(unsigned long iPoint, unsigned short iField) const {
    return connSend[Index[iPoint] + iField];
  }

  /*!
   * \brief Get the Processor ID a Point belongs to.
   * \param[in] iPoint - global renumbered ID of the point
   * \return The rank/processor number.
   */
  virtual unsigned short FindProcessor(unsigned long iPoint) const {
    return linearPartitioner.GetRankContainingIndex(iPoint);
  }

  /*!
   * \brief Get the vector containing the names of the output fields
   * \return Vector of strings containing the field names
   */
  const vector<string>& GetFieldNames() const{
    return fieldNames;
  }

  /*!
   * \brief Get the spatial dimension
   * \return The spatial dimension
   */
  unsigned short GetnDim() const {
    return nDim;
  }

  /*!
   * \brief Set the total number of elements after sorting individual element types
   */
  void SetTotalElements();

};
