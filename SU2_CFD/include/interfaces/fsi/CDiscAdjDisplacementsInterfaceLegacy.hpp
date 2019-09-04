/*!
 * \file CDiscAdjDisplacementsInterfaceLegacy.hpp
 * \brief Declaration and inlines of the class to transfer structural displacements
 *        from a structural zone into a fluid zone in a discrete adjoint simulation.
 * \author Ruben Sanchez
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

#include "../CInterface.hpp"

class CDiscAdjDisplacementsInterfaceLegacy : public CInterface {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjDisplacementsInterfaceLegacy(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjDisplacementsInterfaceLegacy(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CDiscAdjDisplacementsInterfaceLegacy(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                             CGeometry *donor_geometry, CGeometry *target_geometry,
                             CConfig *donor_config, CConfig *target_config);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  void GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
                         unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
                          CConfig *flow_config, unsigned long Marker_Flow,
                          unsigned long Vertex_Flow, unsigned long Point_Flow);

};

