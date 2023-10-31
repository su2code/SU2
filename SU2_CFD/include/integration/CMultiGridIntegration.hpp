/*!
 * \file CMultiGridIntegration.hpp
 * \brief Declaration of class for time integration using a multigrid method.
 * \author F. Palacios, T. Economon
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

#include "CIntegration.hpp"

/*!
 * \class CMultiGridIntegration
 * \ingroup Drivers
 * \brief Class for time integration using a multigrid method.
 * \author F. Palacios
 */
class CMultiGridIntegration final : public CIntegration {
public:
  /*!
   * \brief Constructor of the class.
   */
  CMultiGridIntegration();

  /*!
   * \brief This subroutine calls the MultiGrid_Cycle and also prepare the multigrid levels and the monitoring.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                           CNumerics ******numerics_container, CConfig **config,
                           unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) override;

private:
  /*!
   * \brief Perform a Full-Approximation Storage (FAS) Multigrid.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] mu - Variable for controlling the kind of multigrid algorithm.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void MultiGrid_Cycle(CGeometry ****geometry, CSolver *****solver_container,
                       CNumerics ******numerics_container, CConfig **config,
                       unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
                       unsigned short iZone, unsigned short iInst);

  /*!
   * \brief Compute the forcing term.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                       CGeometry *geo_coarse, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Add the truncation error to the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] flow - Flow solution.
   */
  void SetResidual_Term(CGeometry *geometry, CSolver *flow);

  /*!
   * \brief Set the value of the corrected fine grid solution.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Compute the gradient in coarse grid using the fine grid information.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                              CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Compute the non-dimensional parameters.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   * \param[in] monitor - value of the non-dimensional parameters for monitoring the convergence.
   */
  void NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container,
                                 CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem,
                                 su2double *monitor);

  /*!
   * \brief Compute the fine solution from a coarse solution.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                               CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Compute the fine grid correction from the coarse solution.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                 CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Do an implicit smoothing of the prolongated correction.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] solution - Container vector with all the solutions on the finest grid.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nSmooth - Number of smoothing iterations.
   * \param[in] val_smooth_coeff - Relaxation factor.
   * \param[in] config - Definition of the particular problem.
   */
  void SmoothProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                    unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config);

  /*!
   * \brief Restrict solution from fine grid to a coarse grid.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] InclSharedDomain - Include the shared domain in the interpolation.
   */
  void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                              CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Initialize the adjoint solution using the primal problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void Adjoint_Setup(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                     unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);

};
