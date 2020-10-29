/*!
 * \file CMeshFEM_DG.hpp
 * \brief Class definition for a mesh object for the DG-FEM solver.
 *        The implementations are in the <i>CMeshFEM_DG.cpp</i> file.
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

#include "CMeshFEM_Base.hpp"
#include "CVolumeElementFEM_DG.hpp"

using namespace std;

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains a grid for the DG FEM solver.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CMeshFEM_DG: public CMeshFEM_Base {
protected:
  unsigned long nVolElemTot{0};    /*!< \brief Total number of local volume elements, including halos. */

  vector<unsigned long> nVolElemOwnedPerTimeLevel;    /*!< \brief Number of owned local volume elements
                                                                  per time level. Cumulative storage. */
  vector<unsigned long> nVolElemInternalPerTimeLevel; /*!< \brief Number of internal local volume elements per
                                                                  time level. Internal means that the solution
                                                                  data does not need to be communicated. */
  vector<unsigned long> nVolElemHaloPerTimeLevel;    /*!< \brief Number of local halo volume elements
                                                                 per time level. Cumulative storage. */

  vector<CVolumeElementFEM_DG> volElem;      /*!< \brief Vector of the local volume elements, including halos. */

public:
  /*!
   * \brief Constructor of the class.
   */
  CMeshFEM_DG(void) : CMeshFEM_Base() {}

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
  ~CMeshFEM_DG(void) {}

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
