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
  su2double **mZeros_Aux  = nullptr;         /*!< \brief Submatrix to make zeros and impose Dirichlet boundary conditions. */
  su2double **mId_Aux = nullptr;             /*!< \brief Diagonal submatrix to impose Dirichelt boundary conditions. */

  unsigned short dir;             /*!< \brief If we separate dimensions this tells us in what dimension we currently are. */

  CSysVector<su2double> auxVecInp; /*!< \brief Auxiliar vectors for output and debugging */
  CSysVector<su2double> auxVecRHS; /*!< \brief Auxiliar vectors for output and debugging */
  CSysVector<su2double> auxVecOut; /*!< \brief Auxiliar vectors for output and debugging */

  CVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  std::vector<su2double> deltaP; /*!< \brief The smoothed gradient with respect to the design variables. */


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
   * \brief Main routine for applying the solver
   */
  void ApplyGradientSmoothing(CGeometry *geometry,
                              CSolver *solver,
                              CNumerics **numerics,
                              CConfig *config);

  /*!
   * \brief Assemble the stiffness matrix
   */
  void Compute_StiffMatrix(CGeometry *geometry,
                           CNumerics **numerics,
                           CConfig *config);

  /*!
   * \brief Calculate the RHS of the PDE
   */
  void Compute_Residual(CGeometry *geometry,
                        CSolver *solver,
                        CConfig *config);

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
                    unsigned short val_marker);

  /*!
   * \brief Set Neumann boundary conditions
   */
  void BC_Neumann(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics **numerics,
                  CConfig *config,
                  unsigned short val_marker);

  /*!
   * \brief Call the linear systems solver
   */
  void Solve_Linear_System(CGeometry *geometry,
                           CConfig *config);

  /*!
   * \brief Extract the solution of the linear solver and store it in the sensitivities of the nodes
   */
  void WriteSensitivities(CGeometry *geometry,
                          CSolver *solver,
                          CConfig *config,
                          unsigned long val_marker=0);

  /*!
   * \brief Get the value of the reference coordinate to set on the element structure.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] indexNode - Index of the node.
   * \param[in] iDim - Dimension required.
   */
  su2double Get_ValCoord(CGeometry *geometry,
                         unsigned long indexNode,
                         unsigned short iDim)  {
    return geometry->nodes->GetCoord(indexNode, iDim);
  }

  /*!
   * \brief Main routine to apply the method only on the surface
   */
  void ApplyGradientSmoothingOnSurface(CGeometry *geometry,
                                       CSolver *solver,
                                       CNumerics **numerics,
                                       CConfig *config,
                                       unsigned long val_marker);

  /*!
   * \brief Compute the stiffness matrix of the surface mesh
   */
  void Compute_Surface_StiffMatrix(CGeometry *geometry,
                                   CNumerics **numerics,
                                   CConfig *config,
                                   unsigned long val_marker);

  /*!
   * \brief Compute the RHS of the PDE on the surface mesh
   */
  void Compute_Surface_Residual(CGeometry *geometry,
                                CSolver *solver,
                                CConfig *config,
                                unsigned long val_marker);

  /*!
   * \brief Set Dirichlet boundary conditions for the surface solver
   */
  void BC_Surface_Dirichlet(CGeometry *geometry,
                            CConfig *config,
                            unsigned short val_marker);

  /*!
   * \brief Extract the Coordinates of the element from geometry
   */
  su2activematrix GetElementCoordinates(CGeometry *geometry,
                                        std::vector<unsigned long>& indexNode,
                                        int EL_KIND = 0);

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
   * \brief Write the content of Sensitivity to the sensitivity in the geometry
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void WriteSens2Geometry(CGeometry *geometry,
                          CConfig *config);

  /*!
   * \brief Read the sensitivity in the geometry
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ReadSens2Geometry(CGeometry *geometry,
                         CConfig *config);

  /*!
   * \brief multiply the surface sensitivities with the parametrization Jacobian
   * \param Jacobian of the parameterization
   * \param bool to clarify if we multiply by transposed or not
   */
  void MultiplyParameterJacobian(su2double *Jacobian,
                                 bool transposed);

  /*!
   * \brief write the DV gradient into a file
   */
  void OutputDVGradient(string out_file="delta_p.txt");

  /*!
   * \brief calculate the original DV gradient similar to SU2_DOT_AD
   */
  void CalculateOriginalGradient(CGeometry *geometry,
                                 CConfig *config,
                                 CVolumetricMovement* grid_movement,
                                 su2double *param_jacobi);

  /*!
   * \brief read or write the surface sensitivity into an Eigen vector
   */
  void WriteReadSurfaceSensitivities(CGeometry *geometry,
                                     CConfig *config,
                                     VectorType& x,
                                     bool write);

  /*!
   * \brief Smooth the system by solving each LES in consecutive order
   */
  void SmoothConsecutive(CGeometry *geometry,
                         CSolver *solver,
                         CNumerics **numerics,
                         CConfig *config,
                         su2double *param_jacobi);

  /*!
   * \brief Return the stiffness matrix
   */
  MatrixType GetStiffnessMatrix(CGeometry *geometry,
                                CNumerics **numerics,
                                CConfig *config);

  /*!
   * \brief Return the stiffness matrix
   */
  MatrixType GetSurfaceStiffnessMatrix(CGeometry *geometry,
                                       CNumerics **numerics,
                                       CConfig *config,
                                       unsigned long val_marker);

  /*!
   * \brief Smooth the system by multiplying out the whole system matrix and solving it
   */
  void SmoothCompleteSystem(CGeometry *geometry,
                            CSolver *solver,
                            CNumerics **numerics,
                            CConfig *config,
                            CVolumetricMovement* grid_movement,
                            su2double *param_jacobi);

};
