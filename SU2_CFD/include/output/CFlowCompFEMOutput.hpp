/*!
 * \file output_flow_fem.hpp
 * \brief Headers of the main subroutines for generating the file outputs.
 *        The subroutines and functions are in the <i>output_structure.cpp</i> file.
 * \author F. Palacios, T. Economon, M. Colonno
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

#include "CFlowOutput.hpp"

class CVariable;

/*! \class CFlowFEMOutput
 *  \brief Output class for compressible Flow problems.
 *  \author R. Sanchez, T. Albring.
 *  \date May 30, 2018.
 */
class CFlowCompFEMOutput : public CFlowOutput {
private:
  
  unsigned short nVar;

  unsigned short turb_model;
  
  bool grid_movement;
  
  su2double RefDensity, RefPressure, RefVel2, factor, RefArea;

public:

  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CFlowCompFEMOutput(CConfig *config, unsigned short nDim);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFlowCompFEMOutput(void);

  /*!
   * \brief Set the history file header
   * \param[in] config - Definition of the particular problem.
   */
  void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver);
  
  /*!
   * \brief LoadSurfaceData
   * \param config
   * \param geometry
   * \param solver
   * \param iPoint
   * \param iMarker
   * \param iVertex
   */
  void LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex);  
  
  /*!
   * \brief SetVolumeOutputFields
   * \param config
   */
  void SetVolumeOutputFields(CConfig *config);
  
  /*!
   * \brief LoadVolumeDataFEM
   * \param config
   * \param geometry
   * \param solver
   * \param iElem
   * \param index
   * \param dof
   */
  void LoadVolumeDataFEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iElem, unsigned long index, unsigned short dof);
  
  /*!
   * \brief SetHistoryOutputFields
   * \param config
   */
  void SetHistoryOutputFields(CConfig *config);

  /*!
   * \brief GetQ_Criterion
   * \param config
   * \param geometry
   * \param node_flow
   * \return 
   */
  su2double GetQ_Criterion(CConfig *config, CGeometry *geometry, CVariable *node_flow);
  
  /*!
   * \brief SetInit_Residuals
   * \param config
   * \return 
   */
  bool SetInit_Residuals(CConfig *config);
  
  /*!
   * \brief SetUpdate_Averages
   * \param config
   * \param dualtime
   * \return 
   */
  bool SetUpdate_Averages(CConfig *config);
  
};
