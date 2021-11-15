/*!
 * \file CGradientSmoothingSolver.hpp
 * \brief SOlver class for Sobolev smoothing of sensitivities.
 * \author T. Dick
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

#include "CSolver.hpp"
#include "../../../Common/include/geometry/elements/CElement.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"

/** minimal include of Eigen for linear system */
#include <Eigen/Dense>
#include <Eigen/Sparse>

using MatrixType = Eigen::Matrix<su2double, Eigen::Dynamic, Eigen::Dynamic>;
using VectorType = Eigen::Matrix<su2double, Eigen::Dynamic, 1>;
using QRdecomposition = Eigen::ColPivHouseholderQR<MatrixType>;
// define the output format you want, you only need one instance of this...
const static Eigen::IOFormat CSVFormat(15, 0, ", ", "\n");

/** Introduction of a new data type to allow compilation with forward mode.
  *
  * This is done for compatibility to the treatment of Jacobian and System in CSolver.hpp.
  * Note that the compuations done dere are always 'passive', i.e. not intended to be differentiated. We only need to define functions depending on this once.
  * Move to Common/include/basic_types later if possible.
  */
#ifndef CODI_FORWARD_TYPE
  using su2matvecscalar = su2mixedfloat;
#else
  using su2matvecscalar = su2double;
#endif

/*! \class CGradientSmoothingSolver
 *  \brief Main class for defining a gradient smoothing.
 *  \author T. Dick.
 *  \date March 25, 2019.
 */
class CGradientSmoothingSolver : public CSolver {
public:
  enum : size_t {MAXNNODE = 8};
  enum : size_t {MAXNDIM = 3};

  unsigned long nElement;

  CElement*** element_container  = nullptr;  /*!< \brief Container which stores the element information. */

  su2double **Jacobian_block = nullptr;      /*!< \brief Submatrix to assemble the Jacobian matrix. */
  su2double **mId_Aux = nullptr;             /*!< \brief Diagonal identity matrix to set blocks in the Jacobian. */

  unsigned int dir;                          /*!< \brief If we separate dimensions this tells us in what dimension we currently are. */

  CSysVector<su2double> auxVec;              /*!< \brief Auxiliar vectors for output and debugging */

  CSysVector<su2double> activeCoord;         /*!< \brief Auxiliar vector to keep the indeces of geometry->vertex->Coord */

  CSysVector<su2matvecscalar> helperVecIn;   /*!< \brief Helper vectors for projection and matrix vector product (must be su2mixedfloat) */
  CSysVector<su2matvecscalar> helperVecOut;  /*!< \brief Helper vectors for projection and matrix vector product (must be su2mixedfloat) */
  CSysVector<su2matvecscalar> helperVecAux;     /*!< \brief Helper vectors for matrix vector product if working on surface (smaller dim) */

  CVariable* nodes = nullptr;                /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  std::vector<su2double> deltaP;             /*!< \brief The smoothed gradient with respect to the design variables. */

  MatrixType hessian;                        /*!< \brief The approximated Hessian with respect to the design variables. */

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
                                    CConfig *config) override;

  /*!
   * \brief Main routine to apply the method only on the surface for mesh sensitivities
   *        Projects and smoothes only in the normal direction!
   */
  void ApplyGradientSmoothingSurface(CGeometry *geometry,
                                       CSolver *solver,
                                       CNumerics **numerics,
                                       CConfig *config,
                                       unsigned long val_marker) override;

  /*!
   * \brief All steps required for smoothing the whole system on DV level in an iterative way
   */
  void ApplyGradientSmoothingDV(CGeometry *geometry,
                                CSolver *solver,
                                CNumerics **numerics,
                                CSurfaceMovement *surface_movement,
                                CVolumetricMovement *grid_movement,
                                CConfig *config,
                                vector<su2double> additionalGrad) override;


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
  CSysMatrixVectorProduct<su2matvecscalar> GetStiffnessMatrixVectorProduct(CGeometry *geometry,
                                                                           CNumerics **numerics,
                                                                           CConfig *config);

  /*!
   * \brief calculate the original DV gradient similar to SU2_DOT_AD
   */
  void CalculateOriginalGradient(CGeometry *geometry,
                                 CVolumetricMovement* grid_movement,
                                 CConfig *config);

  /*!
   * \brief Return the original gradient for application in OneShot.
   */
  void RecordTapeAndCalculateOriginalGradient(CGeometry *geometry,
                                              CSurfaceMovement *surface_movement,
                                              CVolumetricMovement *grid_movement,
                                              CConfig *config) override;

  /*!
   * \brief Return current parameter gradient.
   */
  inline vector<su2double> GetDeltaP() override { return deltaP; }

  /*!
   * \brief Return a handle for the Hessian.
   */
  void GetHessianMatrix() override {}

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
                       CSysVector<su2matvecscalar>& result,
                       CSysVector<su2double>& registeredCoord,
                       CConfig *config);

  /*!
   * \brief Reverse evaluate parameterization Jacobian.
   * \param geometry
   * \param config
   * \param surface_movement
   */
  void ProjectMeshToDV(CGeometry *geometry,
                       CSysVector<su2matvecscalar>& sensitivity,
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
                      CConfig *config,
                      CSolver *solver) override;

  /*!
   * \brief Store smoothed sensitivties back into the adjoint solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The solver container holding all terms of the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void OutputSensitivity(CGeometry *geometry,
                         CConfig *config,
                         CSolver *solver) override;

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
                            CSysVector<su2matvecscalar> &vector);

  /*!
   * \brief Copy sensitivities from the geometry into a vector
   */
  void ReadVector2Geometry(CGeometry *geometry,
                           CConfig *config,
                           CSysVector<su2matvecscalar> &vector);

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
