/*!
 * \file CVolumetricMovement.hpp
 * \brief Headers of the CVolumetricMovement class.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "CGridMovement.hpp"

/*!
 * \class CVolumetricMovement
 * \brief Class for moving the volumetric numerical grid.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 */
class CVolumetricMovement : public CGridMovement {
 protected:
  unsigned short nDim; /*!< \brief Number of dimensions. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CVolumetricMovement();

  /*!
   * \brief Constructor of the Class.
   */

  CVolumetricMovement(CGeometry* geometry);

  /*!
   * \brief Destructor of the class.
   */
  ~CVolumetricMovement() override;

  /*!
   * \brief Update the dual grid after the grid movement (edges and control volumes).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateDualGrid(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Update the coarse multigrid levels after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateMultiGrid(CGeometry** geometry, CConfig* config);

  /*!
   * \brief Compute the shape functions for hexahedron
   * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
   */
  su2double GetHexa_Volume(su2double CoordCorners[8][3]) const;

  /*!
   * \brief Compute the shape functions for hexahedron
   * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
   */
  su2double GetTetra_Volume(su2double CoordCorners[8][3]) const;

  /*!
   * \brief Compute the shape functions for hexahedron
   * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
   */
  su2double GetPrism_Volume(su2double CoordCorners[8][3]) const;

  /*!
   * \brief Compute the shape functions for hexahedron
   * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
   */
  su2double GetPyram_Volume(su2double CoordCorners[8][3]) const;

  /*!
   * \brief Compute the shape functions for hexahedron
   * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
   */
  su2double GetTriangle_Area(su2double CoordCorners[8][3]) const;

  /*!
   * \brief Compute the shape functions for hexahedron
   * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
   */
  su2double GetQuadrilateral_Area(su2double CoordCorners[8][3]) const;

  /*!
   * \brief Check for negative volumes (all elements) after performing grid deformation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Screen_Output - determines if text is written to screen
   */
  void ComputeDeforming_Element_Volume(CGeometry* geometry, su2double& MinVolume, su2double& MaxVolume,
                                       bool Screen_Output);

  /*!
   * \brief Compute amount of nonconvex elements
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Screen_Output - determines if text is written to screen
   */
  void ComputenNonconvexElements(CGeometry* geometry, bool Screen_Output);

  /*!
   * \brief Unsteady grid movement using rigid mesh rotation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
   */
  void Rigid_Rotation(CGeometry* geometry, CConfig* config, unsigned short iZone, unsigned long iter);

  /*!
   * \brief Unsteady pitching grid movement using rigid mesh motion.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
   */
  void Rigid_Pitching(CGeometry* geometry, CConfig* config, unsigned short iZone, unsigned long iter);

  /*!
   * \brief Unsteady plunging grid movement using rigid mesh motion.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
   */
  void Rigid_Plunging(CGeometry* geometry, CConfig* config, unsigned short iZone, unsigned long iter);

  /*!
   * \brief Unsteady translational grid movement using rigid mesh motion.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
   */
  void Rigid_Translation(CGeometry* geometry, CConfig* config, unsigned short iZone, unsigned long iter);

  /*!
   * \brief Scale the volume grid by a multiplicative factor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   */
  void SetVolume_Scaling(CGeometry* geometry, CConfig* config, bool UpdateGeo);

  /*!
   * \brief Translate the volume grid by a specified displacement vector.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   */
  void SetVolume_Translation(CGeometry* geometry, CConfig* config, bool UpdateGeo);

  /*!
   * \brief Rotate the volume grid around a specified axis and angle.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   */
  void SetVolume_Rotation(CGeometry* geometry, CConfig* config, bool UpdateGeo);

  /*!
   * \brief Grid deformation using the spring analogy method.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
   */
  inline virtual void SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo,
                                            bool Derivative = false, bool ForwardProjectionDerivative = false){};

  /*!
   * \brief Grid deformation using the spring analogy method.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
   */
  inline virtual void SetVolume_Deformation_Elas(CGeometry* geometry, CConfig* config, bool UpdateGeo,
                                                 bool screen_output, bool Derivative = false) {}

  /*!
   * \brief Set the boundary dependencies in the mesh side of the problem
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Boundary_Dependencies(CGeometry** geometry, CConfig* config) {}

  /*!
   * \brief Returns true if a marker should be deformed.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Index of marker in "marker all".
   */
  inline bool IsDeformationMarker(const CConfig* config, unsigned short iMarker) {
    /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
     deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
     meshes after imposing design variable surface deformations (DV_MARKER). ---*/
    const auto KindSU2 = config->GetKind_SU2();
    const bool cfd = KindSU2 == SU2_COMPONENT::SU2_CFD;
    const bool defDot = KindSU2 == SU2_COMPONENT::SU2_DEF || KindSU2 == SU2_COMPONENT::SU2_DOT;
    const bool dvMarker = config->GetMarker_All_DV(iMarker) == YES;
    const bool moving = config->GetMarker_All_Moving(iMarker) == YES;
    return (cfd && ((config->GetDirectDiff() == D_DESIGN && dvMarker) || moving)) || (defDot && dvMarker);
  }
};
