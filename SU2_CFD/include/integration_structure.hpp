/*!
 * \file integration_structure.hpp
 * \brief Headers of the main subroutines for space and time integration. 
 *        The subroutines and functions are in the <i>integration_structure.cpp</i>, 
 *        <i>integration_time.cpp</i>, and <i>integration_notime.cpp</i> files.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "solver_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \class CIntegration
 * \brief Main class for doing the space integration, time integration, and monitoring 
 *        of a system of Partial Differential Equations (PDE).
 * \author F. Palacios
 */
class CIntegration {
protected:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
  su2double Cauchy_Value,  /*!< \brief Summed value of the convergence indicator. */
  Cauchy_Func;      /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter;  /*!< \brief Number of elements of the Cauchy serial. */
  su2double *Cauchy_Serie;      /*!< \brief Complete Cauchy serial. */
  su2double Old_Func,  /*!< \brief Old value of the objective function (the function which is monitored). */
  New_Func;      /*!< \brief Current value of the objective function (the function which is monitored). */
  bool Convergence,    /*!< \brief To indicate if the flow solver (direct, adjoint, or linearized) has converged or not. */
  Convergence_FSI,    /*!< \brief To indicate if the FSI problem has converged or not. */
  Convergence_FullMG;    /*!< \brief To indicate if the Full Multigrid has converged and it is necessary to add a new level. */
  su2double InitResidual;  /*!< \brief Initial value of the residual to evaluate the convergence level. */

public:
  
  /*! 
   * \brief Constructor of the class. 
   */
  CIntegration(CConfig *config);
  
  /*! 
   * \brief Destructor of the class. 
   */
  virtual ~CIntegration(void);
  
  /*! 
   * \brief Do the space integration of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Space_Integration(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, 
               unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);

  /*!
   * \brief Do the space integration of the numerical system on a FEM framework.
   * \author R. Sanchez
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Space_Integration_FEM(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config,
               unsigned short RunTime_EqSystem);

  /*! 
   * \brief Do the time integration (explicit or implicit) of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
              unsigned short iRKStep, unsigned short RunTime_EqSystem);
  
  /*! 
   * \brief Do the time integration (explicit or implicit) of the numerical system on a FEM framework.
   * \author R. Sanchez
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void Time_Integration_FEM(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config,
                  unsigned short RunTime_EqSystem);

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

  /*! 
   * \brief Average the scalar output in case there is a unsteady solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iteration - Current iteration.
   * \param[in] monitor - Objective function to be averaged.
   * \param[in] iMesh - Index of the mesh.
   */
  void Average_Monitoring(CGeometry *geometry, CConfig *config,
                unsigned long Iteration, su2double monitor, unsigned short iMesh);

  /*!
   * \brief Get the value of the convergence.
   * \return Level of convergence of the solution.
   */
  su2double GetCauchy_Value(void);
  
  /*! 
   * \brief Get the indicator of the convergence for the direct, adjoint and linearized problem.
   * \return <code>TRUE</code> means that the convergence criteria is satisfied; 
   *         otherwise <code>FALSE</code>.
   */
  bool GetConvergence(void);
  
  /*! 
   * \brief Get the indicator of the convergence for the Fluid-Structure Interaction problem.
   * \return <code>TRUE</code> means that the convergence criteria is satisfied;
   *         otherwise <code>FALSE</code>.
   */
  bool GetConvergence_FSI(void);

  /*!
   * \brief Set the indicator of the convergence.
   * \param[in] value - <code>TRUE</code> means that the convergence criteria is satisfied; 
   *            otherwise <code>FALSE</code>.
   */
  void SetConvergence(bool value);
  

  /*!
   * \brief Set the indicator of the convergence for FSI.
   * \param[in] valueFSI - <code>TRUE</code> means that the convergence criteria for FSI is satisfied;
   *            otherwise <code>FALSE</code>.
   */
  void SetConvergence_FSI(bool valueFSI);


  /*! 
   * \brief Get the indicator of the convergence for the full multigrid problem.
   * \return <code>TRUE</code> means that the convergence criteria is satisfied; 
   *         otherwise <code>FALSE</code>.
   */
  bool GetConvergence_FullMG(void);
  
  /*!
   * \brief Save the solution, and volume at different time steps.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetDualTime_Solver(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned short iMesh);
  
  /*! 
   * \brief Save the structural solution at different time steps.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Structural solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetStructural_Solver(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Save the structural solution at different time steps.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Structural solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEM_StructuralSolver(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                  CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);
  
  /*! 
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] mu - Variable for controlling the kind of multigrid algorithm.   
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void MultiGrid_Cycle(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                 CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
                 unsigned short iZone, unsigned short iInst);
  
  /*! 
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] monitor - value of the non-dimensional parameters for monitoring the convergence.
   */
  virtual void NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container, 
                                         CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem,
                                         su2double *monitor);
  
  /*! 
   * \brief A virtual member.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh);

  /*! 
   * \brief A virtual member.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, 
                     CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
  
  /*! 
   * \brief A virtual member.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
                    CGeometry *geo_coarse, CConfig *config);
  
  /*! 
   * \brief A virtual member.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] InclSharedDomain - Include the shared domain in the interpolation.
   */
  virtual void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] InclSharedDomain - Include the shared domain in the interpolation.
   */
  virtual void SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*! 
   * \brief A virtual member.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                    CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
  
  /*! 
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] flow - Flow solution.
   */
  virtual void SetResidual_Term(CGeometry *geometry, CSolver *flow);
  
  /*! 
   * \brief A virtual member.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, 
                 CConfig *config, unsigned short iMesh);
  
  /*! 
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                  CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);


  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                  CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);

  
  /*! 
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  virtual void SetPotential_Solver(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                                   CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone);
  
  /*!
   * \brief A virtual member.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] solution - Container vector with all the solutions on the finest grid.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nSmooth - Number of smoothing iterations.
   * \param[in] val_smooth_coeff - Relaxation factor.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Smooth_Solution(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                       unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config);

};

/*!
 * \class CMultiGridIntegration
 * \brief Class for doing the numerical integration using a multigrid method.
 * \author F. Palacios
 */
class CMultiGridIntegration : public CIntegration {
protected:
    
public:
  
  /*! 
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CMultiGridIntegration(CConfig *config);
  
  /*! 
   * \brief Destructor of the class. 
   */
  ~CMultiGridIntegration(void);
  
  /*! 
   * \brief This subroutine calls the MultiGrid_Cycle and also prepare the multigrid levels and the monitoring.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
               CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);
  
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
  void MultiGrid_Cycle(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                     CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);
  
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
  void GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
                                 CGeometry *geo_coarse, CConfig *config);
  
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
   * \brief Do an implicit smoothing of the solution.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] solution - Container vector with all the solutions on the finest grid.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nSmooth - Number of smoothing iterations.
   * \param[in] val_smooth_coeff - Relaxation factor.
   * \param[in] config - Definition of the particular problem.
   */
  void Smooth_Solution(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                    unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config);

  /*!
   * \brief Set the value of the corrected fine grid solution.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh);

  /*! 
   * \brief Compute truncation error in the coarse grid using the fine grid information. 
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
                CGeometry *geo_coarse, CConfig *config);

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
  void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
  
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
   * \brief Add the truncation error to the residual. 
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] flow - Flow solution.
   */
  void SetResidual_Term(CGeometry *geometry, CSolver *flow);
  
  /*! 
   * \brief Compute the forcing term. 
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, 
             CConfig *config, unsigned short iMesh);
};

/*! 
 * \class CSingleGridIntegration
 * \brief Class for doing the numerical integration of the turbulence model.
 * \author A. Bueno.
 */
class CSingleGridIntegration : public CIntegration {
public:
  
  /*! 
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CSingleGridIntegration(CConfig *config);
  
  /*! 
   * \brief Destructor of the class. 
   */
  ~CSingleGridIntegration(void);
  
  /*! 
   * \brief Do the numerical integration (implicit) of the turbulence solver. 
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
               CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);
  
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
  void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
  
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
  void SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
 
};


/*!
 * \class CStructuralIntegration
 * \brief Class for doing the numerical integration of the structural model.
 * \author R. Sanchez.
 */
class CStructuralIntegration : public CIntegration {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CStructuralIntegration(CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CStructuralIntegration(void);

  /*!
   * \brief Do the numerical integration (implicit) of the structural solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
               CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);

};

/*!
 * \class CFEM_DG_Integration
 * \brief Class for integration with the FEM DG solver.
 * \author E. van der Weide, T. Economon
 * \version 6.2.0 "Falcon"
 */
class CFEM_DG_Integration : public CIntegration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_DG_Integration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DG_Integration(void);
  
  /*!
   * \brief Do the numerical integration (implicit) of the turbulence solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                            CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst);
  /*!
   * \brief Perform the spatial integration of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iStep - Current step of the Runge-Kutta iteration for the RK schemes
                        and the step in the local time stepping for ADER-DG.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Space_Integration(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem);
  
  /*!
   * \brief Perform the time integration (explicit or implicit) of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iStep - Current step of the Runge-Kutta iteration for the RK schemes
                        and the step in the local time stepping for ADER-DG.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iStep, unsigned short RunTime_EqSystem);
};

#include "integration_structure.inl"
