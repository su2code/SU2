#pragma once

#include "CParallelDataSorter.hpp"
#include "CFEMDataSorter.hpp"

class CSurfaceFEMDataSorter : public CParallelDataSorter{
  
  CFEMDataSorter* volume_sorter;
  std::vector<unsigned long> globalSurfaceDOFIDs;
    
public:
  
  CSurfaceFEMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields, CFEMDataSorter* volume_sorter);
  
  ~CSurfaceFEMDataSorter();
  
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
  void SortConnectivity(CConfig *config, CGeometry *geometry,  bool val_sort);
  
  /*!
   * \brief Get the global index of the surface point
   * \param Local surface index
   * \return Global index
   */
  unsigned long GetGlobalIndex(unsigned long iPoint) { return globalSurfaceDOFIDs[iPoint]; }
  
private:
  
  /*!
   * \brief Sort the connectivity for a single surface element type into a linear partitioning across all processors (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

};
