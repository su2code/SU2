/*!
 * \file CMeshFEM_DG.hpp
 * \brief Class definition for a mesh object for the DG-FEM solver.
 *        The implementations are in the <i>CMeshFEM_FG.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CMeshFEM.hpp"

using namespace std;

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains a grid for the DG FEM solver.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CMeshFEM_DG: public CMeshFEM {

public:
  /*!
   * \brief Constructor of the class.
   */
  CMeshFEM_DG(void) : CMeshFEM() {}

  /*!
   * \overload
   * \brief Redistributes the grid over the ranks and creates the halo layer.
   * \param[in] geometry - The linear distributed grid that must be redistributed.
   * \param[in] config   - Definition of the particular problem.
   */
  CMeshFEM_DG(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CMeshFEM_DG(void) {}

  /*!
   * \brief Function to compute the coordinates of the integration points.
   */
  void CoordinatesIntegrationPoints(void);

  /*!
   * \brief Function to compute the coordinates of solution DOFs.
   */
  void CoordinatesSolDOFs(void);

  /*!
   * \brief Function to create the faces used in the DG formulation.
   * \param[in] config - Definition of the particular problem.
   */
  void CreateFaces(CConfig *config);

  /*!
   * \brief Function to create the standard volume elements.
   * \param[in] config - Definition of the particular problem.
   */
  void CreateStandardVolumeElements(CConfig *config);

  /*!
   * \brief Function to compute the grid velocities for static problems.
   * \param[in] config             - Definition of the particular problem.
   * \param[in] Kind_Grid_Movement - The type of prescribed grid motion.
   * \param[in] iZone              - The currently active zone number.
   */
  void InitStaticMeshMovement(CConfig              *config,
                              const unsigned short Kind_Grid_Movement,
                              const unsigned short iZone);

  /*!
   * \brief Function, which computes a length scale of the volume elements.
   * \param[in] config - Definition of the particular problem.
   */
  void LengthScaleVolumeElements(void);

  /*!
   * \brief Function, which computes the metric terms of the surface
   *        elements, both internal faces and physical boundary faces.
   * \param[in] config - Definition of the particular problem.
   */
  void MetricTermsSurfaceElements(CConfig *config);

  /*!
   * \brief Function, which computes the metric terms of the
   *        volume elements.
   * \param[in] config - Definition of the particular problem.
   */
  void MetricTermsVolumeElements(CConfig *config);

  /*!
   * \brief Set the local index that correspond with the global numbering index.
   */
  void SetGlobal_to_Local_Point() override;

  /*!
   * \brief Function, which carries out the preprocessing tasks when wall functions are used.
   * \param[in] config - Definition of the particular problem.
   */
  void WallFunctionPreprocessing(CConfig *config);

protected:
};

/*!
 * \class CDummyMeshFEM_DG
 * \brief Class for defining a DG geometry that does not contain any points/elements.
 *        Can be used for initializing other classes that depend on the geometry without
 *        going through the time-consuming mesh initialization and paritioning.
 * \author T. Albring
 */
class CDummyMeshFEM_DG : public CMeshFEM_DG {

public:
  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CDummyMeshFEM_DG(CConfig *config);

};
