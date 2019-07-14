#pragma once

#include "../../../Common/include/mpi_structure.hpp"
#include "../../../Common/include/option_structure.hpp"

class CGeometry;
class CConfig;

class CParallelDataSorter{
protected:
  
  int rank; 
  int size;
  
  unsigned long nGlobal_Poin_Par;   // Global number of nodes with halos
  unsigned long nGlobal_Elem_Par;  // Global number of elems without halos
  unsigned long nParallel_Poin;
  unsigned long nParallel_Line,
  nParallel_Tria,
  nParallel_Quad,
  nParallel_Tetr,
  nParallel_Hexa,
  nParallel_Pris,
  nParallel_Pyra;
  int *Conn_Line_Par;
  int *Conn_Tria_Par;  // triangle 1 = Conn_Tria[0], Conn_Tria[1], Conn_Tria[3]
  int *Conn_Quad_Par;
  int *Conn_Tetr_Par;
  int *Conn_Hexa_Par;
  int *Conn_Pris_Par;
  int *Conn_Pyra_Par;
  
  unsigned long nGlobalPoint_Sort;
  unsigned long nLocalPoint_Sort;
  unsigned long nPoint_Restart;

  unsigned long *beg_node;
  unsigned long *end_node;

  unsigned long *nPoint_Lin;
  unsigned long *nPoint_Cum;
  
  unsigned short GlobalField_Counter;
  
  su2double** Parallel_Data;
  
  bool connectivity_sorted;  
  
public:
  
  CParallelDataSorter(CConfig *config, unsigned short nFields);
  
  virtual ~CParallelDataSorter();
  
  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  virtual void SortOutputData(CConfig *config, CGeometry *geometry){}

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing.
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
  unsigned long GetNodeBegin(unsigned short rank){return beg_node[rank];}
  
  /*!
   * \brief Ending node ID of the linear partition owned by a specific processor.
   * \input rank - the processor rank.
   * \return The ending node ID.
   */
  unsigned long GetNodeEnd(unsigned short rank){return end_node[rank];}
  
  /*!
   * \brief Get the value of the linear partitioned data.
   * \input iField - the output field ID.
   * \input iPoint - the point ID.
   * \return the value of the data field at a point.
   */
  su2double GetData(unsigned short iField, unsigned long iPoint) {return Parallel_Data[iField][iPoint];}
  
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
  unsigned long GetnPointCumulative(unsigned short rank){return nPoint_Cum[rank];}
  
  /*!
   * \brief Get the linear number of points
   * \input rank - the processor rank.
   * \return The linear number of points up to certain processor rank.
   */
  unsigned long GetnPointLinear(unsigned short rank){return nPoint_Lin[rank];}  
  
  /*!
   * \brief Check whether the current connectivity is sorted (i.e. if SortConnectivity has been called)
   * \return <TRUE> if the connectivity is sorted.
   */  
  bool GetConnectivitySorted(){return connectivity_sorted;}

  /*!
   * \brief Deallocate temporary memory needed for merging and writing output data in parallel.
   */
  void DeallocateData();
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing connectivity in parallel.
   * \param[in] surf_sol - if <TRUE>, surface connectivity is deallocated, otherwise the volume connectivity.
   */
  void DeallocateConnectivity();
  
  /*!
   * \brief CreateLinearPartition
   */
  void CreateLinearPartition(unsigned long nGlobalPoint);
  
  /*!
   * \brief FindProcessor
   * \param global_index
   * \return 
   */
  unsigned short FindProcessor(unsigned long global_index);
  
};
