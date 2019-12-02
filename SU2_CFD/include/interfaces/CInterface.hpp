/*!
 * \file CInterface.hpp
 * \brief Declarations and inlines of the transfer structure.
 *        The subroutines and functions are in the physics folders.
 * \author R. Sanchez
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

#include "../../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "../../../Common/include/config_structure.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../solver_structure.hpp"

using namespace std;

/*!
 * \class CInterface
 * \brief Main class for defining the physical transfer of information.
 * \author R. Sanchez
 * \version 7.0.0 "Blackbird"
 */

class CInterface {
protected:

  int rank,  /*!< \brief MPI Rank. */
  size;      /*!< \brief MPI Size. */

  su2double *Physical_Constants;
  su2double *Donor_Variable;
  su2double *Target_Variable;
  bool valAggregated;

  /*--- Mixing Plane interface variable ---*/
  su2double *SpanValueCoeffTarget;
  unsigned short *SpanLevelDonor;
  unsigned short nSpanMaxAllZones;

  unsigned short nVar;

public:
  /*!
   * \brief Constructor of the class.
   */
  CInterface(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] val_nConst - Number of physical constants that need to be taken into account.
   * \param[in] config - Definition of the particular problem.
   */
  CInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CInterface(void);

  /*!
   * \brief Interpolate data and broadcast it into all processors, for nonmatching meshes.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void BroadcastData(CSolver *donor_solution, CSolver *target_solution,
                     CGeometry *donor_geometry, CGeometry *target_geometry,
                     CConfig *donor_config, CConfig *target_config);
  /*!
   * \brief A virtual member.
   */

  inline virtual void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                                            CGeometry *donor_geometry, CGeometry *target_geometry,
                                            CConfig *donor_config, CConfig *target_config) { }
  /*!
   * \brief A virtual member.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  inline virtual void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                        CConfig *donor_config, unsigned long Marker_Donor,
                                        unsigned long Vertex_Donor, unsigned long Point_Donor) { }

  /*!
   * \brief Initializes the target variable.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] nDonorPoints - Number of donor points.
   */
  inline virtual void InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                                                unsigned long Vertex_Target, unsigned short nDonorPoints){
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Target_Variable[iVar] = 0.0;
  }

  /*!
   * \brief Recovers the target variable from the buffer of su2doubles that was broadcasted.
   * \param[in] indexPoint_iVertex - index of the vertex in the buffer array.
   * \param[in] Buffer_Bcast_Variables - full broadcasted buffer array of doubles.
   * \param[in] donorCoeff - value of the donor coefficient.
   */
  inline virtual void RecoverTarget_Variable(long indexPoint_iVertex, su2double *Buffer_Bcast_Variables,
                                             su2double donorCoeff){
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Target_Variable[iVar] += donorCoeff * Buffer_Bcast_Variables[indexPoint_iVertex*nVar+iVar];

  }

  /*!
   * \brief A virtual member.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  inline virtual void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                         CConfig *target_config, unsigned long Marker_Target,
                                         unsigned long Vertex_Target, unsigned long Point_Target) { }

  /*!
   * \brief A virtual member.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_zone      - Index of the donorZone.
   */

  inline virtual void SetAverageValues(CSolver *donor_solution, CSolver *target_solution,
                                       unsigned short donorZone){ }

  /*!
   * \brief A virtual member.
   * \param[in] donor_geometry - Geometry of the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_zone      - Index of the donorZone.
   */
  inline virtual void SetAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry,
                                               unsigned short donorZone) { }

  /*!
   * \brief A virtual member.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  inline virtual void SetSpanWiseLevels(CConfig *donor_config, CConfig *target_config) { }

  /*!
   * \brief Transfer pre-processing for the mixing plane inteface.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void PreprocessAverage(CGeometry *donor_geometry, CGeometry *target_geometry,
                         CConfig *donor_config, CConfig *target_config, unsigned short iMarkerInt);

  /*!
   * \brief Interpolate data and scatter it into different processors, for matching meshes.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void AllgatherAverage(CSolver *donor_solution, CSolver *target_solution,
                        CGeometry *donor_geometry, CGeometry *target_geometry,
                        CConfig *donor_config, CConfig *target_config, unsigned short iMarkerInt);

  /*!
   * \brief Interpolate data and scatter it into different processors, for matching meshes.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GatherAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone);

  /*!
   * \brief Exchange Average geometrical value beteween zones .
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GatherAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone);

};

