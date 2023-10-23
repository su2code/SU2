/*!
 * \file CFlowTractionInterface.hpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone.
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
#include <unordered_map>

/*!
 * \brief Fluid-structure interface (forces or tractions).
 * \ingroup Interfaces
 */
class CFlowTractionInterface : public CInterface {
protected:
  const bool conservative;
  std::unordered_map<unsigned long, su2double> vertexArea;

  /*!
   * \brief Sets the dimensional factor for pressure and the consistent_interpolation flag.
   * \param[in] flow_config - Definition of the fluid (donor) problem.
   * \param[in] struct_config - Definition of the structural (target) problem.
   * \param[in] geometry - FEA geometry.
   * \param[in] solution - FEA solver.
   */
  void Preprocess(const CConfig *flow_config, const CConfig *struct_config, CGeometry *geometry, CSolver *solution);

  /*!
   * \brief Computes vertex areas (FEA side) for when tractions need to be integrated.
   * \param[in] config - Definition of the structural (target) problem.
   * \param[in] geometry - FEA geometry.
   * \param[in] solution - FEA solver.
   */
  void ComputeVertexAreas(const CConfig *config, CGeometry *geometry, CSolver *solution);

public:
  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] val_nConst - Number of constants.
   * \param[in] config - Definition of the particular problem.
   * \param[in] integrate_tractions_ - Whether to integrate the fluid tractions
   *            (to transfer forces when using conservative interpolation).
   */
  CFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst,
                         const CConfig *config, bool integrate_tractions_);

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
                             const CConfig *donor_config, const CConfig *target_config) override;

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  void GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, const CConfig *flow_config,
                         unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Flow) override;

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                          const CConfig *fea_config, unsigned long Marker_Struct,
                          unsigned long Vertex_Struct, unsigned long Point_Struct) override;

};
