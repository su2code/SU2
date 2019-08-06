#pragma once

#include "CParallelDataSorter.hpp"
#include <vector>

class CFVMDataSorter : public CParallelDataSorter{
  
private:
  std::vector<std::vector<su2double> >* Local_Data;
  int* Local_Halo;
  
public:
  CFVMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields, std::vector<std::vector<su2double> >& Local_Data);
  
  ~CFVMDataSorter();
  
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
  
  unsigned long GetGlobalIndex(unsigned long iPoint) { return beg_node[rank] + iPoint; }
  
  bool GetHalo(unsigned long iPoint){return Local_Halo[iPoint];}
  
private:
  
  void SetHaloPoints(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Sort the connectivity for a single volume element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  void SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type, bool val_sort);
  
};
