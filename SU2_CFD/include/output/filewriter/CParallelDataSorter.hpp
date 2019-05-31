#pragma once
#include "../../../Common/include/geometry_structure.hpp"

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
  map<unsigned long,unsigned long> Global2Renumber, 
                                   Renumber2Global;
  
  unsigned long nGlobalPoint_Sort;
  unsigned long nLocalPoint_Sort;
  unsigned long nPoint_Restart;
  int *Local_Halo_Sort;

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
  virtual void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort){}
 
  unsigned long GetnPoints(){return nParallel_Poin;}
  
  unsigned long GetnPointsGlobal(){return nGlobal_Poin_Par;}
  
  unsigned long GetnElem(){return nGlobal_Elem_Par;}
    
  unsigned long GetnElem(GEO_TYPE type);
    
  unsigned long GetElem_Connectivity(GEO_TYPE type, unsigned long iElem, unsigned long iNode);
  
  unsigned long GetNodeBegin(unsigned short rank){return beg_node[rank];}
  
  unsigned long GetNodeEnd(unsigned short rank){return end_node[rank];}
  
  su2double GetData(unsigned short iField, unsigned long iPoint) {return Parallel_Data[iField][iPoint];}
  
  virtual unsigned long GetGlobalIndex(unsigned long iPoint){return 0;}
  
  unsigned long GetnPointCumulative(){return nPoint_Cum[rank];}
  
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
  
};
