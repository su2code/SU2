/*!
 * \file CInriaFileWriter.hpp
 * \brief Headers fo the GMF file writer class.
 * \author B. Munguía
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include "../../../include/solver_structure.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include "../../../Common/include/config_structure.hpp"

class CInriaFileWriter final: public CFileWriter{

  public:
  
  /*!
   * \brief File extension
   */
  const static string fileExt;
  
  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   * \param[in] fileName - The name of the file
   * \param[in] data_sorter - The parallel sorted data to write
   */  
  CInriaFileWriter(vector<string> fields, unsigned short nDim, 
                 string fileName, CParallelDataSorter* data_sorter);
  
  /*!
   * \brief Destructor
   */
  ~CInriaFileWriter() override;

  /*!
   * \brief Write a native GMF solb file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Flow, adjoint or linearized solution.
   * \param[in] val_iZone - iZone index.
   */
  void WriteInriaRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);

  /*!
   * \brief Write a native GMF meshb file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void WriteInriaMesh(CConfig *config, CGeometry *geometry);

  // /*!
  //  * \brief Organizes the all the solutions for AMG.
  //  * \param[in] solver_container - Container vector with all the solutions.
  //  * \param[in] geometry - Geometrical definition of the problem.
  //  * \param[in] config - Definition of the particular problem.
  //  * \param[in] val_nZone - Total number of domains in the grid file.
  //  */
  // void SetResult_Parallel(CSolver *****solver_container, CGeometry ****geometry, CConfig **config,
  //                         unsigned short val_nZone);

  // /*!
  //  * \brief Returns a sorted solution value for AMG for the current processor.
  //  */
  // vector<vector<passivedouble> > GetResult_Parallel(void);

  // /*!
  //  * \brief Cleans up the sorted solutions.
  //  */
  // void CleanResult_Parallel(void);

  // /*!
  //  * \brief Organizes the connectivity for AMG.
  //  * \param[in] geometry - Geometrical definition of the problem.
  //  * \param[in] config - Definition of the particular problem.
  //  * \param[in] val_nZone - Total number of domains in the grid file.
  //  */
   
  // void SetConnectivity_Parallel(CGeometry ****geometry, CConfig **config, unsigned short val_nZone);

  // /*!
  //  * \brief Get the connectivity of edges.
  //  */
  // vector<vector<unsigned long> > GetConnEdg(CConfig *config, CGeometry *geometry);

  // /*!
  //  * \brief Get the connectivity of tris.
  //  */
  // vector<vector<unsigned long> > GetConnTri(CConfig *config, CGeometry *geometry);

  // /*!
  //  * \brief Get the connectivity of tets.
  //  */
  // vector<vector<unsigned long> > GetConnTet(CConfig *config, CGeometry *geometry);

  // /*!
  //  * \brief Cleans up the connectivity.
  //  */
  // void CleanConnectivity_Parallel(void);

};