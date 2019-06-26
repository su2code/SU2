#pragma once

#include "CParallelDataSorter.hpp"

class CFEMDataSorter : public CParallelDataSorter{
private:
  
  std::vector<std::vector<su2double> >* Local_Data;
  
  
public:
  CFEMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields, std::vector<std::vector<su2double> > &Local_Data);
  
  ~CFEMDataSorter();
  
  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing connectivity in parallel.
   * \param[in] surf_sol - if <TRUE>, surface connectivity is deallocated, otherwise the volume connectivity.
   */
  void DeallocateConnectivity_Parallel(bool surf_sol);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing output data in parallel.
   */
  void DeallocateData_Parallel();
  
  unsigned long GetGlobalIndex(unsigned long iPoint) { return beg_node[rank] + iPoint; }
  
private:
  
  
  /*!
   * \brief Sort the connectivity for a single volume element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);
  
};
