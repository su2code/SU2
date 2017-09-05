/*!
 * \file hybrid_anisotropy_model.hpp
 * \brief 
 * \author C. Pederson
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#include "stdio.h"
#include "math.h"

// Forward declarations to resolve circular dependencies
class CSolver;

using namespace std;

/*!
 * \class CHybrid_Visc_Anisotropy
 * \brief Base (abstract) class for the subgrid anisotropy model for the
 *        turbulent stress.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Visc_Anisotropy {
 protected:
  su2double** eddy_visc_anisotropy;
  const unsigned short nDim;
 public:

  /**
   * \brief Default constructor for base class
   * \param[in] nDim - Number of dimensions of the problem (e.g. 2D or 3D)
   */
  CHybrid_Visc_Anisotropy(unsigned short nDim);

  /**
   * \brief Base class destructor
   */
  virtual ~CHybrid_Visc_Anisotropy();

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  virtual void CalculateViscAnisotropy() = 0;

  /**
   * \brief Retrieves the turbulent stress anisotropy tensor.
   * \return The turbulent stress anisotropy.
   */
  su2double** GetViscAnisotropy();

  /**
   * \brief Sets a rank-2 tensor used in the anisotropy model.
   *
   * In the SGS anisotropy, this is the approximate structure function. This is
   * left abstract, in order to allow different possible tensors for different
   * implementations.
   *
   * \param[in] val_tensor - The tensor used in the model
   */
  virtual void SetTensor(su2double** val_tensor) = 0;

  /**
   * \brief Sets a scalar used in the anisotropy model.
   *
   * In the SGS anisotropy, this is the resolution inadequacy and/or
   * the RANS weighting parameter.  This is left abstract, in order to allow
   * different possible scalars for different implementations.
   *
   * \param val_scalars - The scalars for the model
   */
  virtual void SetScalars(vector<su2double> val_scalars) = 0;
};

/*!
 * \class CHybrid_Isotropic_Visc
 * \brief Subgrid anisotropy model based on the approximate 2nd order
 *        structure function.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Isotropic_Visc : public CHybrid_Visc_Anisotropy {
 public:

  CHybrid_Isotropic_Visc(unsigned short nDim);

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  void CalculateViscAnisotropy();
  
  /**
   * \brief This method does nothing. No tensors are needed.
   * \param val_tensor - This value is not used.
   */
  void SetTensor(su2double** val_tensor) {}

  /**
   * \brief This method does nothing. No scalars are needed.
   * \param val_scalar
   */
  void SetScalars(vector<su2double> val_scalars) {}
};

/*!
 * \class CHybrid_Aniso_Q
 * \brief Subgrid anisotropy model based on the approximate 2nd order
 *        structure function.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Aniso_Q : public CHybrid_Visc_Anisotropy {
 protected:
  su2double** Qstar; ///< \brief The approximate two-point structure function at the grid resolution
  su2double   rans_weight; ///< \brief The weight given to the isotropic (RANS) portion

 public:
  /**
   * \brief Constructor for the anisotropy model.
   * \param[in] nDim - The number of dimensions (2D or 3D)
   */
  CHybrid_Aniso_Q(unsigned short nDim);
  
  /**
   * \brief Sets the approximate structure function.
   * \param[in] val_approx_struct_func
   */
  void SetTensor(su2double** val_approx_struct_func);

  /**
   * \brief Sets the RANS weight
   * \param[in] val_scalars
   */
  void SetScalars(vector<su2double> val_scalars);

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  void CalculateViscAnisotropy();
};



/*!
 * \class CAbstract_Hybrid_Mediator
 * \brief Base abstract class for a hybrid model mediator object.
 *
 * In order to decouple the RANS model, the subgrid model, the hybrid parameter,
 * and the resolved flow, a mediator object is necessary.  This allows the
 * RANS, subgrid, hybrid parameter, and resolved flow equations to follow the
 * single responsibility principle, while this class makes sure they
 * have the information they need.
 *
 * The main purpose of this class is to pass variables to where they
 * need to go.
 *
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CAbstract_Hybrid_Mediator {
 public:
  /*!
   * \brief A virtual destructor.
   */
  virtual ~CAbstract_Hybrid_Mediator() {};

  /**
   * \brief Retrieve and pass along all necessary info for the RANS model.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] rans_numerics - The source numerics for the turb. solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  virtual void SetupRANSNumerics(CGeometry* geometry,
                                 CSolver **solver_container,
                                 CNumerics* rans_numerics,
                                 unsigned short iPoint,
                                 unsigned short jPoint) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for the hybrid solver
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The node being evaluated
   */
  virtual void SetupHybridParamSolver(CGeometry* geometry,
                                      CSolver **solver_container,
                                      unsigned short iPoint) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for the hybrid parameter numerics
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] hybrid_numerics - The source numerics for the hybrid parameter
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  virtual void SetupHybridParamNumerics(CGeometry* geometry,
                                        CSolver **solver_container,
                                        CNumerics *hybrid_numerics,
                                        unsigned short iPoint,
                                        unsigned short jPoint) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for the stress
   *        anisotropy model.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] hybrid_anisotropy - The hybrid anisotropy model
   * \param[in] iPoint - The node being evaluated
   */
  virtual void SetupStressAnisotropy(CGeometry* geometry,
                                     CSolver **solver_container,
                                     CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                                     unsigned short iPoint) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for resolved numerics.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] visc_numerics - The viscous numerics for the resolved flow solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  virtual void SetupResolvedFlowNumerics(CGeometry* geometry,
                                     CSolver **solver_container,
                                     CNumerics* visc_numerics,
                                     unsigned short iPoint,
                                     unsigned short jPoint) = 0;
};


/*!
 * \class CHybrid_Mediator
 * \brief Mediator object for the Q-based hybrid RANS/LES model.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Mediator : public CAbstract_Hybrid_Mediator {
 protected:

  unsigned short nDim;
  const su2double C_sf; /*!> \brief Model constant relating the structure function to the unresolved turbulent kinetic energy  */
  su2double C_zeta; /*!> \brief Scaling constant for the transformation tensor zeta */
  su2double **Q,        /*!> \brief An approximate 2nd order structure function tensor */
            **Qapprox;  /*!> \brief An approximate 2nd order structure function tensor (used for temporary calculations) */
  std::vector<std::vector<su2double> > constants;
  CConfig* config;

  /*--- Data structures for LAPACK ---*/
#ifdef HAVE_LAPACK
  int info, lwork;
  su2double wkopt;
  su2double* work;
  su2double eigval[3], eigvec[9], vr[9], wi[3];
  su2double mat[9];
  int num_found;
  int isupp[3];
#endif

  /*!
   * \brief Calculates the resolution inadequacy parameter
   * \param[in] Q - The approximate 2nd order structure function
   * \param[in] v2 - The v2 value from Durbin's k-eps-v2-f model
   * \return The resolution inadequacy parameter
   */
  su2double CalculateRk(su2double** Q, su2double v2);

  /*!
   * \brief Projects the resolution on a specific vector
   * \param[in] resolution_tensor - The tensor representing separation distances
   * \param[in] direction - The direction vector (assumed to be normalized)
   * \return The magnitude of the resolution tensor projected on the direction.
   */
  su2double GetProjResolution(su2double** resolution_tensor,
                              vector<su2double> direction);

  /**
   * \brief Uses a resolution tensor and a gradient-gradient tensor to build an
   *        approximate two-point structure function tensor
   * \param[in] val_ResolutionTensor - A tensor representing cell-cell distances
   * \param[in] val_PrimVar_Grad - The gradient in the resolved velocity field.
   * \param[out] val_Q - An approximate resolution-scale two-point second-order
   *                 structure function.
   */
  template <class T>
  void CalculateApproxStructFunc(T val_ResolutionTensor,
                                 su2double** val_PrimVar_Grad,
                                 su2double** val_Q);

  /*!
   * \brief Loads the model fit constants from *.dat files
   * \param[in] filename - The base name for the files (e.g. [filename]0.dat)
   * \return The 3 sets of constants pulled from the files.
   */
  vector<vector<su2double> > LoadConstants(string filename);

  /*!
   * \brief Solve for the eigenvalues of Q, given eigenvalues of M
   *
   * Solves for the eigenvalues of the expected value of the contracted velocity
   * differences at grid resolution (Q), given the eigenvalues of the
   * resolution tensor.
   *
   * \param[in] eig_values_M - Eigenvalues of the resolution tensor.
   * \return The eigenvalues of the expected value of the approx SF tensor.
   */
  vector<su2double> GetEigValues_Q(vector<su2double> eig_values_M);

  /*!
   * \brief Calculates the eigenvalues of a modified resolution tensor.
   * \param eig_values_M - Eigenvalues of the grid-based resolution tensor.
   * \return Eigenvalues of a transformation mapping the approximate velocity
   *         differences at grid resolution to the two-point second-order
   *         structure function.
   */
  vector<su2double> GetEigValues_Zeta(vector<su2double> eig_values_M);

 public:
  /*!
   * \brief Builds a transformation for the approximate structure function.
   * \param[in] values_M - The cell-to-cell distances in the "principal
   *            "directions"
   */
  vector<vector<su2double> > BuildZeta(su2double* values_M,
                                       su2double** vectors_M);


  /**
   * \brief Constructor for the hybrid mediator object.
   * \param[in] nDim - The number of dimensions of the problem
   * \param[in] CConfig - The configuration for the current zone
   */
  CHybrid_Mediator(int nDim, CConfig* config, string filename="");

  /**
   * \brief Destructor for the hybrid mediator object.
   */
  ~CHybrid_Mediator();

  /**
   * \brief RANS needs the hybrid parameter (the energy flow parameter).
   *
   * This function sets the hybrid parameter from the previous timestep.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] rans_numerics - The source numerics for the turb. solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupRANSNumerics(CGeometry* geometry, CSolver **solver_container,
                         CNumerics* rans_numerics,
                         unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The hybrid solver needs the resolution adequacy parameter, which
   *        is dependent on RANS results.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The node being evaluated
   */
  void SetupHybridParamSolver(CGeometry* geometry, CSolver **solver_container,
                           unsigned short iPoint);

  /**
   * \brief The hybrid numerics need the turbulence length and timescales as
   *        well as the resolution adequacy parameter.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] hybrid_numerics - The source numerics for the hybrid parameter
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupHybridParamNumerics(CGeometry* geometry, CSolver **solver_container,
                                CNumerics *hybrid_param_numerics,
                                unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The turbulent stress anisotropy needs the weighting factor and the
   *        second order structure function.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] hybrid_anisotropy - The hybrid anisotropy model
   * \param[in] iPoint - The node being evaluated
   */
  void SetupStressAnisotropy(CGeometry* geometry,
                             CSolver **solver_container,
                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                             unsigned short iPoint);

  /**
   * \brief The resolved flow needs the hybrid parameter (the energy flow
   *        parameter) and the turbulent stress anisotropy tensor.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] visc_numerics - The viscous numerics for the resolved flow solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupResolvedFlowNumerics(CGeometry* geometry,
                             CSolver **solver_container,
                             CNumerics* visc_numerics,
                             unsigned short iPoint,
                             unsigned short jPoint);

  /**
   * \brief Returns the constants for the numerical fit for the resolution tensor.
   * \return Constants for the numerical fit for the resolution tensor.
   */
  vector<vector<su2double> > GetConstants();


  void SolveEigen(su2double** M, vector<su2double> &eigvalues,
                  vector<vector<su2double> > &eigvectors);
};

/*!
 * \class CHybrid_Dummy_Mediator
 * \brief Mediator object for RANS-only operation; isotropic stress and dummy hybrid parameter
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Dummy_Mediator : public CAbstract_Hybrid_Mediator {
 protected:

  unsigned short nDim;
  su2double* dummy_alpha; // A default value of alpha to pass around

 public:

  /**
   * \brief Constructor for the hybrid mediator object.
   * \param[in] nDim - The number of dimensions of the problem
   * \param[in] CConfig - The configuration for the current zone
   */
  CHybrid_Dummy_Mediator(int nDim, CConfig* config);

  /**
   * \brief Destructor for the hybrid mediator object.
   */
  ~CHybrid_Dummy_Mediator();

  /**
   * \brief RANS needs the hybrid parameter (the energy flow parameter).
   *
   * This function sets the hybrid parameter from the previous timestep.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] rans_numerics - The source numerics for the turb. solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupRANSNumerics(CGeometry* geometry, CSolver **solver_container,
                         CNumerics* rans_numerics,
                         unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The hybrid solver doesn't need anything.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The node being evaluated
   */
  void SetupHybridParamSolver(CGeometry* geometry, CSolver **solver_container,
                           unsigned short iPoint);

  /**
   * \brief The hybrid numerics don't need anything.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] hybrid_numerics - The source numerics for the hybrid parameter
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupHybridParamNumerics(CGeometry* geometry, CSolver **solver_container,
                                CNumerics *hybrid_param_numerics,
                                unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The turbulent stress anisotropy needs the weighting factor and the
   *        second order structure function.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] hybrid_anisotropy - The hybrid anisotropy model
   * \param[in] iPoint - The node being evaluated
   */
  void SetupStressAnisotropy(CGeometry* geometry,
                             CSolver **solver_container,
                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                             unsigned short iPoint);

  /**
   * \brief The resolved flow needs the hybrid parameter (the energy flow
   *        parameter) and the turbulent stress anisotropy tensor.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] visc_numerics - The viscous numerics for the resolved flow solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupResolvedFlowNumerics(CGeometry* geometry,
                             CSolver **solver_container,
                             CNumerics* visc_numerics,
                             unsigned short iPoint,
                             unsigned short jPoint);
};

/*--- Template definitions:
 * These must be placed with the template declarations.  They can be kept here
 * or moved to an *.inl file that is included in this header. ---*/

template <class T>
void CHybrid_Mediator::CalculateApproxStructFunc(T val_ResolutionTensor,
                                                 su2double** val_PrimVar_Grad,
                                                 su2double** val_Q) {
  unsigned int iDim, jDim, kDim, lDim, mDim;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      val_Q[iDim][jDim] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      for (kDim = 0; kDim < nDim; kDim++)
        for (lDim = 0; lDim < nDim; lDim++)
          for (mDim = 0; mDim < nDim; mDim++)
            val_Q[iDim][jDim] += val_ResolutionTensor[iDim][mDim]*
                             val_PrimVar_Grad[kDim+1][mDim]*
                             val_PrimVar_Grad[kDim+1][lDim]*
                             val_ResolutionTensor[lDim][jDim];
}

