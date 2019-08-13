#pragma once

#include "CParallelDataSorter.hpp"
#include "CFVMDataSorter.hpp"

class CSurfaceFVMDataSorter : public CParallelDataSorter{
  
  CFVMDataSorter* volume_sorter;
  map<unsigned long,unsigned long> Renumber2Global;
  
public:
  CSurfaceFVMDataSorter(CConfig *config, CGeometry* geometry, unsigned short nFields, CFVMDataSorter* volume_sorter);
  
  ~CSurfaceFVMDataSorter();
  
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
   * \param[in] surf - boolean controlling whether surface <TRUE> or volume connectivity <FALSE> should be sorted.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort);
  
  /*!
   * \brief Get the global index of the surface point
   * \param Local surface index
   * \return Global index
   */  
  unsigned long GetGlobalIndex(unsigned long iPoint) { return Renumber2Global[iPoint]; }
  
  
private:
  
  /*!
   * \brief Sort the connectivity for a single surface element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

};
