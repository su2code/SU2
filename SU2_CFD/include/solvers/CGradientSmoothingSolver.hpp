/*!
 * \file CGradientSmoothingSolver.hpp
 * \brief SOlver class for Sobolev smoothing of sensitivities.
 * \author T. Dick
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

#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"
#include "../variables/CSobolevSmoothingVariable.hpp"
#include "CFEASolverBase.hpp"

/*!
 * \brief Main class for defining a Sobolev-based gradient smoothing.
 * \author T. Dick.
 * \ingroup GradSmooth
 */
class CGradientSmoothingSolver final : public CFEASolverBase {
public:

/** Introduction of a new alias for the data type to allow compilation with forward mode.
 *
 * This is done for compatibility to the treatment of Jacobian and System in CSolver.hpp.
 * Note that the computations done here are always 'passive', i.e. not intended to be differentiated.
 * We only need to define functions depending on this once.
 */
#ifndef CODI_FORWARD_TYPE
  typedef su2mixedfloat su2matvecscalar;
#else
  typedef su2double su2matvecscalar;
#endif

 private:
  unsigned int curDim;                       /*!< \brief If we separate dimensions this tells us in what dimension we currently are. */

  CSysVector<su2double> activeCoord;         /*!< \brief Auxiliar vector to keep the indeces of geometry->vertex->Coord */

  CSysVector<su2matvecscalar> helperVecIn;   /*!< \brief Helper vectors for projection and matrix vector product (must be su2mixedfloat) */
  CSysVector<su2matvecscalar> helperVecOut;  /*!< \brief Helper vectors for projection and matrix vector product (must be su2mixedfloat) */
  CSysVector<su2matvecscalar> helperVecAux;  /*!< \brief Helper vectors for matrix vector product if working on surface (smaller dim) */

  std::vector<su2double> deltaP;             /*!< \brief The smoothed gradient with respect to the design variables. */

  CSquareMatrixCM hessian;                   /*!< \brief The approximated Hessian with respect to the design variables. */

  std::vector<bool> visited;                 /*! <\brief Stores already visited points for surface applications with multiple markers. */

  /*!
   * \brief The highest level in the variable hierarchy all derived solvers can safely use,
   * CVariable is the common denominator between the FEA and Mesh deformationd variables.
   */
  CSobolevSmoothingVariable* nodes = nullptr;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

 public:
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
  void ApplyGradientSmoothingVolume(CGeometry* geometry, CNumerics* numerics, const CConfig* config) override;

  /*!
   * \brief Main routine to apply the method only on the surface for mesh sensitivities
   *        Projects and smoothes only in the normal direction!
   */
  void ApplyGradientSmoothingSurface(CGeometry* geometry, CNumerics* numerics, const CConfig* config) override;

  /*!
   * \brief All steps required for smoothing the whole system on DV level in an iterative way
   */
  void ApplyGradientSmoothingDV(CGeometry *geometry,
                                CNumerics *numerics,
                                CSurfaceMovement *surface_movement,
                                CVolumetricMovement *grid_movement,
                                CConfig *config,
                                su2double** Gradient) override;

 private:

  /*!
   * \brief Assemble the stiffness matrix
   */
  void Compute_StiffMatrix(CGeometry* geometry, CNumerics* numerics, const CConfig* config);

  /*!
   * \brief Compute the stiffness matrix of the surface mesh
   */
  void Compute_Surface_StiffMatrix(CGeometry* geometry, CNumerics* numerics, const CConfig* config,
                                   unsigned long val_marker, unsigned int nSurfDim = 1);

  /*!
   * \brief Calculate the RHS of the PDE
   */
  void Compute_Residual(CGeometry* geometry, const CConfig* config);

  /*!
   * \brief Compute the RHS of the PDE on the surface mesh
   */
  void Compute_Surface_Residual(CGeometry* geometry, const CConfig* config, unsigned long val_marker);

  /*!
   * \brief Set the boundary conditions
   */
  void Impose_BC(const CGeometry* geometry, const CConfig* config);

  /*!
   * \brief Set Dirichlet boundary conditions
   */
  void BC_Dirichlet(const CGeometry* geometry, const CConfig* config, unsigned int val_marker);

  /*!
   * \brief Set Dirichlet boundary conditions for the surface solver
   */
  void BC_Surface_Dirichlet(const CGeometry* geometry, const CConfig* config, unsigned int val_marker);

  /*!
   * \brief Call the linear systems solver
   */
  void Solve_Linear_System(CGeometry* geometry, const CConfig* config);

  /*!
   * \brief Get the matrix vector product with the StiffnessMatrix
   * \note This always applies the stiffness matrix for all dimensions independent of each other!
   */
  template <typename scalar_type>
  CSysMatrixVectorProduct<scalar_type> GetStiffnessMatrixVectorProduct(CGeometry* geometry, CNumerics* numerics,
                                                                       const CConfig* config);

  /*!
   * \brief calculate the original DV gradient similar to SU2_DOT_AD
   */
  void CalculateOriginalGradient(CGeometry *geometry,
                                 CVolumetricMovement* grid_movement,
                                 CConfig *config,
                                 su2double **Gradient);

  /*!
   * \brief Return the original gradient for application in OneShot.
   */
  void RecordTapeAndCalculateOriginalGradient(CGeometry *geometry,
                                              CSurfaceMovement *surface_movement,
                                              CVolumetricMovement *grid_movement,
                                              CConfig *config,
                                              su2double** Gradient) override;

  /*!
   * \brief write the DV gradient into a file
   */
  void OutputDVGradient(const CConfig* config, string out_file = "delta_p.txt");

  /*!
   * \brief Record a tape containing the parameter Jacobian.
   * \param geometry - Geometrical definition of the problem.
   * \param surface_movement - Surface movement of the problem.
   * \param registeredCoord - Indexes of the affected coordinates.
   * \param config - Definition of the particular problem.
   */
  void RecordParameterizationJacobian(CGeometry *geometry,
                                      CSurfaceMovement *surface_movement,
                                      CSysVector<su2double> &registeredCoord,
                                      CConfig *config);

  /*!
   * \brief Forward evaluate parameterization Jacobian.
   */
  void ProjectDVtoMesh(CGeometry *geometry,
                       std::vector<su2double>& seeding,
                       CSysVector<su2matvecscalar>& result,
                       CSysVector<su2double>& registeredCoord,
                       CConfig *config);

  /*!
   * \brief Reverse evaluate parameterization Jacobian.
   */
  void ProjectMeshToDV(CGeometry *geometry,
                       CSysVector<su2matvecscalar>& sensitivity,
                       std::vector<su2double>& output,
                       CSysVector<su2double> &registeredCoord,
                       CConfig *config);

  /*!
   * \brief Write the content of sensitivity in the nodes to the sensitivity in the geometry
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void WriteSensToGeometry(CGeometry* geometry) const override;

  /*!
   * \brief Read the sensitivity for the geometry into the nodes
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void ReadSensFromGeometry(const CGeometry *geometry) override;

  /*!
   * \brief Write the solution of the linear solver into the sensitivities of the nodes
   */
  void WriteSensitivity(CGeometry* geometry, const CConfig* config);

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
   * \brief Extra entries to eliminate in the linear system
   */
  void Set_VertexEliminationSchedule(CGeometry* geometry, const CConfig* config);

  /*!
   * \brief Complete the calculation of the surface stiffness matrix
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void Complete_Surface_StiffMatrix(const CGeometry* geometry);

  /*!
   * \brief Get the element container index and number of nodes of a given VTK type.
   * \param[in] VTK_Type - Type of element.
   * \param[out] EL_KIND - Element container index.
   * \param[out] nNodes - Number of nodes.
   */
  static void GetElemKindAndNumNodes(unsigned short VTK_Type, int& EL_KIND, unsigned int& nNodes) {
    /*--- if we need higher order quadrature rules overide some of the element kinds ---*/
    switch (VTK_Type) {
      case LINE:          nNodes = 2; EL_KIND = EL_LINE;  break;
      case TRIANGLE:      nNodes = 3; EL_KIND = EL_TRIA2;  break;
      case QUADRILATERAL: nNodes = 4; EL_KIND = EL_QUAD;  break;
      case TETRAHEDRON:   nNodes = 4; EL_KIND = EL_TETRA2; break;
      case PYRAMID:       nNodes = 5; EL_KIND = EL_PYRAM2; break;
      case PRISM:         nNodes = 6; EL_KIND = EL_PRISM; break;
      case HEXAHEDRON:    nNodes = 8; EL_KIND = EL_HEXA;  break;
      default: assert(false); nNodes = 0; EL_KIND = -(1<<30); break;
    }
  }

  /*!
   * \brief Set the current working dimension, if the seperate dimension option is set.
   * \param[in] iDim - the dimension we are currently working in.
   */
  inline void SetCurrentDim(unsigned int iDim) {
    curDim=iDim;
  }

  /*!
   * \brief Return the current working dimension.
   */
  inline unsigned int GetCurrentDim() const {
    return curDim;
  }

};
