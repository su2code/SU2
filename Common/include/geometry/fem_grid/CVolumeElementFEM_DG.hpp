/*!
 * \file CVolumeElementFEM_DG.hpp
 * \brief Class for a volume element for the DG-FEM solver.
 *        The implementations are in the <i>CVolumeElementFEM_DG.cpp</i> file.
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

#include "CVolumeElementFEM_Base.hpp"
#include "../../toolboxes/classes_multiple_integers.hpp"

using namespace std;

/*!
 * \class CVolumeElementFEM_DG
 * \brief Class to store a volume element for the DG-FEM solver.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CVolumeElementFEM_DG final: public CVolumeElementFEM_Base {
public:

  bool elemAdjLowerTimeLevel = false;  /*!< \brief Whether or not the element is adjacent to a lower
                                                   time level when time accurate local time stepping
                                                   is employed. */

  short periodIndexToDonor = -1;    /*!< \brief The index of the periodic transformation to the donor
                                                element. Only for halo elements. A -1 indicates no
                                                periodic transformation. */
  unsigned short timeLevel = 0;    /*!< \brief Time level of the element when time accurate local
                                               time stepping is employed. */

  unsigned int factTimeLevel = 0;  /*!< \brief Number of local time steps for this element
                                               compared to the largest time step when time
                                               accurate local time stepping is employed. */

  unsigned long offsetDOFsSolGlobal = 0; /*!< \brief Global offset of the solution DOFs of this element. */

  su2double deltaTime = -1.0;        /*!< \brief Time step for this element. */

  vector<unsigned long> internalFaceIDs; /*!< \brief Vector, which contains the ID's of the internal
                                                     faces of this element. */

  vector<CUnsignedLong2T> boundaryFaceIDs; /*!< \brief Vector, which contains the ID's of the boundary
                                                       faces of this element. First index is the boundary,
                                                       second index is the face. */

  ColMajorMatrix<su2double> metricTerms2ndDer; /*!< \brief The metric terms needed for the computation
                                                           of the 2nd derivatives in the integration
                                                           points. Only determined when needed (ADER-DG
                                                           with non-aliased predictor for the
                                                           Navier-Stokes equations). */

  ColMajorMatrix<su2double> solDOFs;                        /*!< \brief The solution in the DOFs. */
  ColMajorMatrix<su2double> solDOFsWork;                    /*!< \brief The working variables in the DOFs. */
  ColMajorMatrix<su2double> solDOFsAux;                     /*!< \brief Auxiliary solution in the DOFs. For the
                                                                        classical RK4 scheme this is the new solution
                                                                        while for ADER it may contain the solution
                                                                        of a previous time level. */
  vector<ColMajorMatrix<su2double> > solDOFsADERPredictor;  /*!< \brief The vector containing the predictor solution
                                                                        in the DOFs (space and time) for ADER. */

  ColMajorMatrix<su2double> dUdVInt;   /*!< \brief The transformation matrix between conservative and
                                                   entropy variables in the volume integration points. ---*/
  
  ColMajorMatrix<su2double> resDOFs;         /*!< \brief The residual in the solution DOFs. */
  ColMajorMatrix<su2double> resTotDOFsADER;  /*!< \brief The total residuals in the solution DOFs for ADER. */

  CFEMStandardElementBase *standardElemFlow = nullptr; /*!< \brief Pointer to the standard element for the
                                                                   standard flow solution variables. */
  CFEMStandardElementBase *standardElemP    = nullptr; /*!< \brief Pointer to the standard element for the
                                                                   pressure for an incompressible flow. */

  /*!
   * \brief Function, which converts the modal solution to a nodal solution
   *        for the flow variables.
   */
  inline void ModalToNodalFlow(void) {standardElemFlow->ModalToNodal(solDOFs);}

  /*!
   * \brief Function, which converts the nodal solution to a modal solution
   *        for the flow variables.
   */
  inline void NodalToModalFlow(void) {standardElemFlow->NodalToModal(solDOFs);}

  /*!
   * \brief Function, which allocate the memory for the compressible
   *        flow variables.
   * \param[in] config - Definition of the particular problem.
   * \param[in] nVar   - Number of flow variables.
   */
  void AllocateCompressibleFlowVar(CConfig        *config,
                                   unsigned short nVar);

  /*!
   * \brief Function, which allocate the memory for the residuals.
   * \param[in] config - Definition of the particular problem.
   * \param[in] nVar   - Number of flow variables.
   */
  void AllocateResiduals(CConfig        *config,
                         unsigned short nVar);

  /*!
   * \brief Function, which computes the gradients of the solution in the integration points.
   * \return  A reference to the gradients of the solution in the integration points.
   */
  vector<ColMajorMatrix<su2double> > &ComputeGradSolIntPoints(void);

  /*!
   * \brief Function, which computes the solution in the integration points.
   * \return  A reference to the solution in the integration points.
   */
  ColMajorMatrix<su2double> &ComputeSolIntPoints(void);

  /*!
   * \brief Function, which adds to the residual the contribution coming
   *        from the multiplication with the basis functions.
   * \param[in] scalarDataInt - The scalar data in the integration points
   *                            to be multiplied by the basis functions.
   */
  void ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt);

  /*!
   * \brief Function, which adds to the residual the contribution coming
   *        from the multiplication with the divergence of the basis functions.
   * \param[in] vectorDataInt - The vector data in the integration points to be
   *                            multiplied by the gradient of the basis functions.
   */
  void ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt);

  /*!
   * \brief Function, which sets the solution in this element
   *        to a constant solution.
   * \param[in] sol      - Solution to be set for this element.
   * \param[in] nVar     - Number of variables to be set.
   * \param[in] startInd - Start index in the solution variable.
   */
  void SetConstantSolution(const su2double *sol,
                           unsigned short  nVar,
                           unsigned short  startInd);
};
