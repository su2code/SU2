/*!
 * \file CMeshFEM_DG.hpp
 * \brief Class definition for a mesh object for the DG-FEM solver.
 *        The implementations are in the <i>CMeshFEM_DG.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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
#include "CInternalFaceFEM_DG.hpp"
#include "../../toolboxes/CSquareMatrixCM.hpp"


/*--- Forward declarations. ---*/
class CFaceOfElement;

using namespace std;

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains a grid for the DG FEM solver.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

  vector<vector<unsigned long> > ownedElemAdjLowTimeLevel; /*!< \brief List of owned elements per time level that are
                                                                       adjacent to elements of the lower time level. */
  vector<vector<unsigned long> > haloElemAdjLowTimeLevel; /*!< \brief List of halo elements per time level that are
                                                                      adjacent to elements of the lower time level. */

  vector<unsigned long> nMatchingFacesInternal;      /*!< \brief Number of matching faces between between two owned elements
                                                                 per time level. Cumulative storage format. */
  vector<unsigned long> nMatchingFacesWithHaloElem;  /*!< \brief Number of matching faces between an owned element and a halo
                                                                 element per time level. Cumulative storage format. */

  vector<CVolumeElementFEM_DG> volElem;       /*!< \brief Vector of the local volume elements, including halos. */
  vector<CInternalFaceFEM_DG>  matchingFaces; /*!< \brief Vector of the local matching internal faces. */

  vector<CFEMStandardElementBase *> standardVolumeElementsSolution;    /*!< \brief Vector of standard volume
                                                                                   elements for the solution. */
  vector<CFEMStandardElementBase *> standardSurfaceElementsSolution;   /*!< \brief Vector of standard surface
                                                                                   elements for the solution. */

  vector<CFEMStandardInternalFaceGrid *> standardInternalFaceGrid;     /*!< \brief Vector of standard elements of internal
                                                                                   matching faces for the grid. */
  vector<CFEMStandardInternalFaceSol *>  standardInternalFaceSolution; /*!< \brief Vector of standard elements of internal
                                                                                   matching faces for the solution. */ 

  map<unsigned long, unsigned long> Global_to_Local_Point; /*!< \brief Global-local mapping for the DOFs. */

  CSquareMatrixCM timeCoefADER_DG;                                      /*!< \brief The time coefficients in the iteration matrix
                                                                                    of the ADER-DG predictor step. */
  ColMajorMatrix<passivedouble> timeInterpolDOFToIntegrationADER_DG;    /*!< \brief The interpolation matrix between the time DOFs
                                                                                    and the time integration points for ADER-DG. */
  ColMajorMatrix<passivedouble> timeInterpolAdjDOFToIntegrationADER_DG; /*!< \brief The interpolation matrix between the time DOFs
                                                                                    of adjacent elements of a higher time level and
                                                                                    the time integration points for ADER-DG. */

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
  ~CMeshFEM_DG(void);

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
   * \brief Set the send receive boundaries of the grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSendReceive(const CConfig *config) override;

  /*!
   * \brief Function, which carries out the preprocessing tasks when wall functions are used.
   * \param[in] config - Definition of the particular problem.
   */
  void WallFunctionPreprocessing(CConfig *config);

private:

  /*!
   * \brief Function, which creates the standard elements for the faces.
   * \param[in] config     - Definition of the particular problem.
   * \param[in] localFaces - Vector, which stores all the faces in the grid.
   */
  void CreateStandardFaces(CConfig                      *config,
                           const vector<CFaceOfElement> &localFaces);

  /*!
   * \brief Function, which creates the standard elements for the solution.
   * \param[in] elemTypes   - Information about the element types to be created.
   * \param[in] locGridDOFs - Location of the grid DOFs, either LGL or equidistant.
   */
  void CreateStandardVolumeElementsSolution(const vector<CUnsignedShort4T> &elemTypes,
                                            const unsigned short           locGridDOFs);

  /*!
   * \brief Set wall distances a specific value
   * \param[in] val - Value to which the wall distance must be set.
   */
  void SetWallDistance(su2double val) override;

  /*!
   * \brief Set the wall distance based on an previously constructed ADT
   * \param[in] config - Definition of the particular problem.
   * \param[in] WallADT - The ADT to compute the wall distance
   */
  void SetWallDistance(const CConfig *config, CADTElemClass* WallADT) override;

  /*!
   * \brief Function, which computes the time coefficients for the ADER-DG predictor step.
   * \param[in] config - Definition of the particular problem.
   */
  void TimeCoefficientsPredictorADER_DG(CConfig *config);
};
