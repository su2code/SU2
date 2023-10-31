/*!
 * \file CConservativeVarsInterface.hpp
 * \brief Declaration and inlines of the class to transfer conservative variables
 *        from a generic zone into another one.
 * \author Ruben Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \brief Transfers conservative variables.
 * \ingroup Interfaces
 */
class CConservativeVarsInterface : public CInterface {
public:
  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   */
  CConservativeVarsInterface(unsigned short val_nVar, unsigned short val_nConst);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   * \param[in] Point_Donor - Index of the donor point.
   */
  void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, const CConfig *donor_config,
                         unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor) override;

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry, const CConfig *target_config,
                          unsigned long Marker_Target, unsigned long Vertex_Target, unsigned long Point_Target) override;


};

