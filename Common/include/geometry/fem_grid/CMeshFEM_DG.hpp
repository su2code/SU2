/*!
 * \file CMeshFEM_DG.hpp
 * \brief Class definition for a mesh object for the DG-FEM solver.
 *        The implementations are in the <i>CMeshFEM_DG.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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
 * \version 7.1.1 "Blackbird"
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

  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

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

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which makes available the number of owned volume elements per time level.
   * \return  The pointer to the data of nVolElemOwnedPerTimeLevel.
   */
  inline unsigned long* GetNVolElemOwnedPerTimeLevel(void) {return nVolElemOwnedPerTimeLevel.data();}

  /*!
   * \brief Function, which makes available the number of internal volume elements per time level.
   * \return  The pointer to the data of nVolElemInternalPerTimeLevel.
   */
  inline unsigned long* GetNVolElemInternalPerTimeLevel(void) {return nVolElemInternalPerTimeLevel.data();}

  /*!
   * \brief Function, which makes available the number of halo volume elements per time level.
   * \return  The pointer to the data of nVolElemHaloPerTimeLevel.
   */
  inline unsigned long* GetNVolElemHaloPerTimeLevel(void) {return nVolElemHaloPerTimeLevel.data();}

  /*!
   * \brief Function, which makes available the vector of vectors containing the owned element
   *        IDs adjacent to elements of a lower time level. Note that a copy is made.
   * \return Copy of ownedElemAdjLowTimeLevel.
   */
  inline vector<vector<unsigned long> > GetOwnedElemAdjLowTimeLevel(void) {return ownedElemAdjLowTimeLevel;}

  /*!
   * \brief Function, which makes available the vector of vectors containing the halo element
   *        IDs adjacent to elements of a lower time level. Note that a copy is made.
   * \return Copy of haloElemAdjLowTimeLevel.
   */
  inline vector<vector<unsigned long> > GetHaloElemAdjLowTimeLevel(void) {return haloElemAdjLowTimeLevel;}

  /*!
  * \brief Function, which makes available the number of matching internal faces
  *        between an owned element and a halo element per time level.
  * \return  The number of matching internal faces between these elements per time level.
  */
  inline unsigned long *GetNMatchingFacesWithHaloElem(void) {return nMatchingFacesWithHaloElem.data();}

 /*!
  * \brief Function, which makes available the number of matching internal faces
  *        between two owned elements per time level.
  * \return  The number of matching internal faces per time level.
  */
  inline unsigned long *GetNMatchingFacesInternal(void) {return nMatchingFacesInternal.data();}

 /*!
  * \brief Function, which makes available the matching internal faces.
  * \return  Pointer to the matching internal faces.
  */
  inline CInternalFaceFEM_DG *GetMatchingFaces(void) {return matchingFaces.data();}

  /*!
   * \brief Function, which makes available the total number of volume elements in the local FEM mesh.
   * \return  Total number of volume elements of the local FEM mesh.
   */
  inline unsigned long GetNVolElemTot(void) const {return nVolElemTot;}

  /*!
   * \brief Function, which makes available the volume elements in the local FEM mesh.
   * \return  Pointer to the volume elements of the local FEM mesh.
   */
  inline CVolumeElementFEM_DG *GetVolElem(void) {return volElem.data();}

  /*!
  * \brief Function, which makes available the time coefficients in the iteration matrix
  *        of the ADER-DG predictor step. Note that a copy is made.
  * \return Copy of the time coefficients in the iteration matrix of ADER-DG.
  */
  inline CSquareMatrixCM GetTimeCoefADER_DG(void) {return timeCoefADER_DG;}

  /*!
   * \brief Function, which makes available the time interpolation matrix between
   *        the time DOFs and time integration points for ADER-DG. Note that a
   *        copy is made.
   * \return Copy of the time interpolation matrix for ADER-DG.
   */
  inline ColMajorMatrix<passivedouble> GetTimeInterpolDOFToIntegrationADER_DG(void) {
    return timeInterpolDOFToIntegrationADER_DG;
  }

  /*!
   * \brief Function, which makes available the time interpolation matrix between the
   *        adjacent time DOFs of the next time level and the time integration points
   *        for ADER-DG. Note that a copy is made.
   * \return Copy of the time interpolation matrix of adjacent time DOFs for ADER-DG.
   */
  inline ColMajorMatrix<passivedouble> GetTimeInterpolAdjDOFToIntegrationADER_DG(void) {
    return timeInterpolAdjDOFToIntegrationADER_DG;
  }

  /*-----------------------------------------------------------------------------------*/
  /*---                      Public member functions.                               ---*/
  /*-----------------------------------------------------------------------------------*/

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
   * \brief Function, which determines the maximum number of DOFs
   *        used in the locally stored grid.
   * \return The maximum number of DOFs.
   */
  unsigned short DetermineMaxNDOFs(void);

  /*!
   * \brief Function, which determines the maximum number of integration
   *        points used in the locally stored grid.
   * \return The maximum number of integration points.
   */
  unsigned short DetermineMaxNIntegration(void);

  /*!
   * \brief Get the local index that correspond with the global numbering index.
   * \param[in] val_ipoint - Global point.
   * \return Local index that correspond with the global index, -1 if not found on the current rank.
   */
  long GetGlobal_to_Local_Point(unsigned long val_ipoint) const override;

  /*!
   * \brief Function to compute the grid velocities for static problems.
   * \param[in] config             - Definition of the particular problem.
   * \param[in] Kind_Grid_Movement - The type of prescribed grid motion.
   * \param[in] iZone              - The currently active zone number.
   */
  void InitStaticMeshMovement(const CConfig        *config,
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

  /*-----------------------------------------------------------------------------------*/
  /*---                      Private member functions.                              ---*/
  /*-----------------------------------------------------------------------------------*/

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
   * \param[in] nAllocVar   - Number of variables for which memory must be allocated
   *                          in the working vectors.
   * \param[in] locGridDOFs - Location of the grid DOFs, either LGL or equidistant.
   */
  void CreateStandardVolumeElementsSolution(const vector<CUnsignedShort4T> &elemTypes,
                                            const unsigned short           nAllocVar,
                                            const unsigned short           locGridDOFs);

  /*!
   * \brief Function, which computes the parametric coordinates of the given
            Cartesian coordinates inside the given parent element.
   * \param[in]  coor           - Cartesian coordinates for which the parametric
                                  coordinates must be determined.
   * \param[in]  parElem        - The high order parent element which contains
                                  the point.
   * \param[in]  subElem        - Low order sub element inside the parent element
                                  which contains the point.
   * \param[in]  weightsSubElem - Interpolation weights inside subElem for the
                                  coordinates. Used for an initial guess.
   * \param[out] parCoor        - Parametric coordinates inside the high order
                                  parent element for the given coordinates.
                                  These parametric coordinates must be computed.
   */
  void HighOrderContainmentSearch(const su2double      *coor,
                                  const unsigned long  parElem,
                                  const unsigned short subElem,
                                  const su2double      *weightsSubElem,
                                  su2double            *parCoor);

  /*!
   * \brief Set wall distances a specific value
   * \param[in] val - Value to which the wall distance must be set.
   */
  void SetWallDistance(su2double val) override;

  /*!
   * \brief Set the wall distance based on an previously constructed ADT
   * \param[in] WallADT - The ADT to compute the wall distance
   * \param[in] config  - Definition of the particular problem.
   * \param[in] iZone   - Zone whose markers made the ADT.
   */
  void SetWallDistance(CADTElemClass* WallADT,
                       const CConfig* config,
                       unsigned short iZone) override;

  /*!
   * \brief Function, which computes the time coefficients for the ADER-DG predictor step.
   * \param[in] config - Definition of the particular problem.
   */
  void TimeCoefficientsPredictorADER_DG(CConfig *config);
};
