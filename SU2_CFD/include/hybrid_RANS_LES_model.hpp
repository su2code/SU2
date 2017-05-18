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
 * \class CHybrid_SGS_Anisotropy
 * \brief Base (abstract) class for the subgrid anisotropy model for the
 *        turbulent stress.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_SGS_Anisotropy{
 protected:
  su2double** stress_anisotropy_tensor;
  const unsigned short nDim;
 public:

  /**
   * \brief Default constructor for base class
   * \param[in] nDim - Number of dimensions of the problem (e.g. 2D or 3D)
   */
  CHybrid_SGS_Anisotropy(unsigned short nDim);

  /**
   * \brief Base class destructor
   */
  ~CHybrid_SGS_Anisotropy();

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  virtual void CalculateStressAnisotropy() = 0;

  /**
   * \brief Retrieves the turbulent stress anisotropy tensor.
   * \return The turbulent stress anisotropy.
   */
  su2double** GetStressAnisotropyTensor();

  // FIXME: These defy the Liskov substitution principle

  virtual void SetApproxStructFunc(su2double** val_approx_struct_func) = 0;
  virtual void SetResolutionAdequacy(su2double val_r_k) = 0;
};

/*!
 * \class CHybrid_Isotropic_Stress
 * \brief Subgrid anisotropy model based on the approximate 2nd order
 *        structure function.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Isotropic_Stress : public CHybrid_SGS_Anisotropy {
 public:

  CHybrid_Isotropic_Stress(unsigned short nDim);

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  void CalculateStressAnisotropy();
  
  void SetApproxStructFunc(su2double** val_approx_struct_func) {}
  void SetResolutionAdequacy(su2double val_r_k) {}
};

/*!
 * \class CHybrid_Aniso_Q
 * \brief Subgrid anisotropy model based on the approximate 2nd order
 *        structure function.
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Aniso_Q : public CHybrid_SGS_Anisotropy {
 protected:
  su2double** Qstar;
  su2double   Qstar_norm;
  su2double   resolution_adequacy;
  su2double CalculateIsotropyWeight(su2double r_k);
 public:
  CHybrid_Aniso_Q(unsigned short nDim);
  
  void SetApproxStructFunc(su2double** val_approx_struct_func);
  void SetResolutionAdequacy(su2double val_r_k);

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  void CalculateStressAnisotropy();

  /**
   * \brief Retrieves the turbulent stress anisotropy tensor.
   * \return The turbulent stress anisotropy.
   */
  su2double** GetStressAnisotropyTensor();
};



/*!
 * \class CAbstract_Hybrid_Mediator
 * \brief Base (abstract) class for a hybrid model mediator object.
 *
 * In order to decouple the RANS model, the subgrid model, the blending,
 * and the mean flow, a mediator object is necessary.  This allows the
 * RANS, subgrid, blending, and mean flow equations to follow the
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
  virtual ~CAbstract_Hybrid_Mediator();

  /**
   * \brief Retrieve and pass along all necessary info for the RANS model.
   */
  virtual void SetupRANSNumerics(CGeometry* geometry,
                                 CSolver **solver_container,
                                 CNumerics* rans_numerics,
                                 unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief Retrieve and pass along all necessary info for the blending.
   */
  virtual void SetupBlendingSolver(CGeometry* geometry,
                                   CSolver **solver_container,
                                   unsigned short iPoint) = 0;

  /**
   *
   * \param geometry
   * \param solver_container
   * \param blending_numerics
   * \param iPoint
   * \param jPoint
   */
  void SetupBlendingNumerics(CGeometry* geometry, CSolver **solver_container,
                             CNumerics *blending_numerics,
                             unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief Retrieve and pass along all necessary info for the stress
   *        anisotropy model.
   */
  virtual void SetupStressAnisotropy(CGeometry* geometry,
                                     CSolver **solver_container,
                                     unsigned short iPoint) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for the resolved flow.
   */
  virtual void SetupMeanFlow(CGeometry* geometry,
                             CSolver **solver_container,
                             CNumerics* visc_numerics,
                             unsigned short iPoint,
                             unsigned short jPoint) = 0;
};


/*!
 * \class CHybrid_Mediator
 * \brief Mediator object for the Q-based hybrid RANS/LES model.
 *
 * In order to decouple the RANS model, the subgrid model, the blending,
 * and the resolved flow, a mediator object is necessary.  This allows the
 * RANS, subgrid, blending, and resolved flow equations to follow the
 * single responsibility principle, while this class makes sure they
 * have the information they need.
 *
 * The main purpose of this class is to pass variables to where they
 * need to go.
 *
 * \author: C. Pederson
 * \version 5.0.0 "Raven"
 */
class CHybrid_Mediator {
 protected:
  /*--- Pointers to the connected numerical objects ---*/
  CHybrid_Aniso_Q* hybrid_anisotropy;

  unsigned short nDim;
  const su2double C_sf; /*!> \brief Smagorinksy constant */
  su2double **Q,
  **ResolutionTensor,
  **PrimVar_Grad_i,
  **PrimVar_Grad_j;
  su2double r_k;

  /*!
   * \brief Calculates the resolution inadequacy parameter
   * \param Q - The approximate 2nd order structure function
   * \param v2 - The v2 value from Durbin's k-eps-v2-f model
   * @return The resolution inadequacy parameter
   */
  su2double CalculateRk(su2double** Q, su2double v2);

  /*!
   * \brief Calculates the approximate 2nd order structure function tensor, Qij
   * \param[in] ResolutionTensor - The resolution (metric) tensor
   * \param[in] PrimVar_Grad - Gradients of the primitive flow variables
   * \param[out] Q - The approximate 2nd order structure function tensor
   */
  void CalculateApproxStructFunc(su2double** ResolutionTensor,
                                 su2double** PrimVar_Grad,
                                 su2double** Q);

 public:

  /**
   * \brief Constructor for the hybrid mediator object.
   * \param[in] nDim - The number of dimensions of the problem
   */
  CHybrid_Mediator(int nDim, CHybrid_Aniso_Q* hybrid_anisotropy);

  /**
   * \brief Destructor for the hybrid mediator object.
   */
  ~CHybrid_Mediator();

  /**
   * Retrieve and calculate all terms which do not depend on any specific model.
   */
  void SetupSelf();

  /**
   * \brief RANS needs the blending coefficient (the energy flow parameter).
   *
   * This function sets the blending coefficient from the previous timestep.
   */
  void SetupRANSNumerics(CGeometry* geometry, CSolver **solver_container,
                         CNumerics* rans_numerics,
                         unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The blending solver needs the resolution adequacy parameter, which
   *        is dependent on RANS results.
   * \param geometry - A pointer to the geometry
   * \param solver_container - An array of solvers
   * \param iPoint - The node being evaluated
   */
  void SetupBlendingSolver(CGeometry* geometry, CSolver **solver_container,
                           unsigned short iPoint);

  /**
   * \brief The blending numerics need the turbulence length and timescales as
   *        well as the resolution adequacy parameter.
   * \param geometry - A pointer to the geometry
   * \param solver_container - An array of solvers
   * \param iPoint - The node being evaluated
   */
  void SetupBlendingNumerics(CGeometry* geometry, CSolver **solver_container,
                             CNumerics *blending_numerics,
                             unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The turbulent stress anisotropy needs the weighting factor and the
   *        second order structure function.
   * \param geometry - A pointer to the geometry
   * \param solver_container - An array of solvers
   * \param iPoint - The node being evaluated
   */
  void SetupStressAnisotropy(CGeometry* geometry,
                             CSolver **solver_container,
                             CHybrid_SGS_Anisotropy* hybrid_anisotropy,
                             unsigned short iPoint);

  /**
   * \brief The resolved flow needs the blending coefficient (the energy flow
   *        parameter) and the turbulent stress anisotropy tensor.
   * \param geometry - A pointer to the geometry
   * \param solver_container - An array of solvers
   * \param iPoint - The node being evaluated
   */
  void SetupMeanFlow(CGeometry* geometry,
                     CSolver **solver_container,
                     CNumerics* visc_numerics,
                     unsigned short iPoint,
                     unsigned short jPoint);
};
