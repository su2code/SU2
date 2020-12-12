/*!
 * \file CGradientSmoothingSolver.hpp
 * \brief SOlver class for Sobolev smoothing of sensitivities.
 * \author T. Dick
 * \version 7.0.5 "Blackbird"
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

#include "CSolver.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"

/*! \class CGradientSmoothingSolver
 *  \brief Main class for defining a gradient smoothing.
 *  \author T. Dick.
 *  \date March 25, 2019.
 */
class CGradientSmoothingSolver : public CSolver {
public:

  unsigned long nElement;

  CElement*** element_container  = nullptr;  /*!< \brief Container which stores the element information. */

  su2double **Jacobian_block = nullptr;      /*!< \brief Submatrix to assemble the Jacobian matrix. */
  su2double **mId_Aux = nullptr;             /*!< \brief Diagonal identity matrix to set blocks in the Jacobian. */

  unsigned int dir;                          /*!< \brief If we separate dimensions this tells us in what dimension we currently are. */

  CSysVector<su2double> auxVec;              /*!< \brief Auxiliar vectors for output and debugging */

  CSysVector<su2double> activeCoord;         /*!< \brief Auxiliar vector to keep the indeces of geometry->vertex->Coord */

  #ifndef CODI_FORWARD_TYPE
    CSysVector<su2mixedfloat> helperVecIn;   /*!< \brief Helper vectors for projection and matrix vector product (must be su2mixedfloat) */
    CSysVector<su2mixedfloat> helperVecOut;  /*!< \brief Helper vectors for projection and matrix vector product (must be su2mixedfloat) */
    CSysVector<su2mixedfloat> helperVecAux;     /*!< \brief Helper vectors for matrix vector product if working on surface (smaller dim) */
  #else
    CSysVector<su2double> helperVecIn;
    CSysVector<su2double> helperVecOut;
    CSysVector<su2double> helperVecAux;     /*!< \brief Helper vectors for matrix vector product if working on surface (smaller dim) */
  #endif

  CVariable* nodes = nullptr;                /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  std::vector<su2double> deltaP;             /*!< \brief The smoothed gradient with respect to the design variables. */

  /*--- Extra vertices for row/column elimination, see Set_VertexEliminationSchedule. ---*/
  vector<unsigned long> ExtraVerticesToEliminate;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

  /*!
   * \brief Constructor of the class.
   */
  CGradientSmoothingSolver(CGeometry *geometry,
                           CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGradientSmoothingSolver(void);

  /*!
   * \brief Main routine for applying the solver on the volume sensitivities
   */
  void ApplyGradientSmoothingVolume(CGeometry *geometry,
                                    CSolver *solver,
                                    CNumerics **numerics,
                                    CConfig *config);

  /*!
   * \brief Main routine to apply the method only on the surface for mesh sensitivities
   *        Projects and smoothes only in the normal direction!
   */
  void ApplyGradientSmoothingSurface(CGeometry *geometry,
                                       CSolver *solver,
                                       CNumerics **numerics,
                                       CConfig *config,
                                       unsigned long val_marker);

  /*!
   * \brief All steps required for smoothing the whole system on DV level in an iterative way
   */
  void ApplyGradientSmoothingDV(CGeometry *geometry,
                                CSolver *solver,
                                CNumerics **numerics,
                                CConfig *config,
                                CSurfaceMovement *surface_movement,
                                CVolumetricMovement *grid_movement);


  /*!
   * \brief Assemble the stiffness matrix
   */
  void Compute_StiffMatrix(CGeometry *geometry,
                           CNumerics **numerics,
                           CConfig *config);

  /*!
   * \brief Compute the stiffness matrix of the surface mesh
   */
  void Compute_Surface_StiffMatrix(CGeometry *geometry,
                                   CNumerics **numerics,
                                   CConfig *config,
                                   unsigned long val_marker,
                                   unsigned int nSurfDim=1);

  /*!
   * \brief Calculate the RHS of the PDE
   */
  void Compute_Residual(CGeometry *geometry,
                        CSolver *solver,
                        CConfig *config);

  /*!
   * \brief Compute the RHS of the PDE on the surface mesh
   */
  void Compute_Surface_Residual(CGeometry *geometry,
                                CSolver *solver,
                                CConfig *config,
                                unsigned long val_marker);

  /*!
   * \brief Set the boundary conditions
   */
  void Impose_BC(CGeometry *geometry,
                 CNumerics **numerics,
                 CConfig *config);

  /*!
   * \brief Set Dirichlet boundary conditions
   */
  void BC_Dirichlet(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics **numerics,
                    CConfig *config,
                    unsigned int val_marker);

  /*!
   * \brief Set Dirichlet boundary conditions for the surface solver
   */
  void BC_Surface_Dirichlet(CGeometry *geometry,
                            CConfig *config,
                            unsigned int val_marker);

  /*!
   * \brief Call the linear systems solver
   */
  void Solve_Linear_System(CGeometry *geometry,
                           CConfig *config);

  /*!
   * \brief Get the matrix vector product with the StiffnessMatrix
   * \note This always applies the stiffness matrix for all dimensions independent of each other!
   */
  CSysMatrixVectorProduct<su2mixedfloat> GetStiffnessMatrixVectorProduct(CGeometry *geometry,                                                                        CNumerics **numerics,
                                                                         CConfig *config);

  /*!
   * \brief calculate the original DV gradient similar to SU2_DOT_AD
   */
  void CalculateOriginalGradient(CGeometry *geometry,
                                 CVolumetricMovement* grid_movement,
                                 CConfig *config);

  /*!
   * \brief write the DV gradient into a file
   */
  void OutputDVGradient(string out_file="delta_p.txt");

  /*!
   * \brief Record a tape containing the parameter Jacobian.
   * \param geometry
   * \param config
   * \param surface_movement
   */
  void RecordParameterizationJacobian(CGeometry *geometry,
                                      CSurfaceMovement *surface_movement,
                                      CSysVector<su2double> &registeredCoord,
                                      CConfig *config);

  /*!
   * \brief Forward evaluate parameterization Jacobian.
   * \param geometry
   * \param config
   * \param surface_movement
   */
  void ProjectDVtoMesh(CGeometry *geometry,
                       std::vector<su2double>& seeding,
                       CSysVector<su2mixedfloat>& result,
                       CSysVector<su2double>& registeredCoord,
                       CConfig *config);

  /*!
   * \brief Reverse evaluate parameterization Jacobian.
   * \param geometry
   * \param config
   * \param surface_movement
   */
  void ProjectMeshToDV(CGeometry *geometry,
                       CSysVector<su2mixedfloat>& sensitivity,
                       std::vector<su2double>& output,
                       CSysVector<su2double> &registeredCoord,
                       CConfig *config);

  /*!
   * \brief Extract and set the sensitivity from the discrete adjoint solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The solver container holding all terms of the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CGeometry *geometry,
                      CSolver **solver,
                      CConfig *config);

  /*!
   * \brief Store smoothed sensitivties back into the adjoint solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The solver container holding all terms of the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void OutputSensitivity(CGeometry *geometry,
                         CSolver **solver,
                         CConfig *config);

  /*!
   * \brief Write the solution of the linear solver into the sensitivities of the nodes
   */
  void WriteSensitivity(CGeometry *geometry,
                          CSolver *solver,
                          CConfig *config,
                          unsigned long val_marker=0);

  /*!
   * \brief Write the content of sensitivity in the nodes to the sensitivity in the geometry
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void WriteSens2Geometry(CGeometry *geometry,
                          CConfig *config);

  /*!
   * \brief Read the sensitivity for the geometry into the nodes
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ReadSens2Geometry(CGeometry *geometry,
                         CConfig *config);

  /*!
   * \brief Copy sensitivities from a vector into the geometry
   */
  void WriteVector2Geometry(CGeometry *geometry,
                            CConfig *config,
                            CSysVector<su2mixedfloat> &vector);

  /*!
   * \brief Copy sensitivities from the geometry into a vector
   */
  void ReadVector2Geometry(CGeometry *geometry,
                           CConfig *config,
                           CSysVector<su2mixedfloat> &vector);

  /*!
   * \brief Get the value of the reference coordinate to set on the element structure.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] indexNode - Index of the node.
   * \param[in] iDim - Dimension required.
   */
  inline su2double Get_ValCoord(CGeometry *geometry,
                         unsigned long indexNode,
                         unsigned int iDim)  {
    return geometry->nodes->GetCoord(indexNode, iDim);
  }

  /*!
   * \brief Extract the Coordinates of the element from geometry
   */
  su2activematrix GetElementCoordinates(CGeometry *geometry,
                                        std::vector<unsigned long>& indexNode,
                                        int EL_KIND = 0);

  /*!
   * \brief Extra entries to eliminate in the linear system
   */
  void Set_VertexEliminationSchedule(CGeometry *geometry,
                                     CConfig *config);

};
