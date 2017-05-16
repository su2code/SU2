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

#include "../../Common/include/mpi_structure.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#include "stdio.h"
#include "math.h"

#include "../SU2_CFD/include/numerics_structure.hpp"

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
   * \brief a virtual destructor
   */
  virtual ~CHybrid_SGS_Anisotropy();

  /**
   * \brief Tells the hybrid model to calculate the turbulent stress anisotropy.
   */
  virtual void CalculateStressAnisotropy() = 0;

  /**
   * \brief Retrieves the turbulent stress anisotropy tensor.
   * @return The turbulent stress anisotropy.
   */
  virtual su2double** GetStressAnisotropyTensor();
};

su2double** CHybrid_SGS_Anisotropy::GetStressAnisotropyTensor() {
  return stress_anisotropy_tensor;
}

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
};

CHybrid_Isotropic_Stress::CHybrid_Isotropic_Stress(unsigned short nDim)
: nDim(nDim) {
  stress_anisotropy_tensor = new su2double*[nDim];
  for (unsigned short iDim=0; iDim < nDim; iDim++) {
    stress_anisotropy_tensor[iDim] = new su2double[nDim];
    for (unsigned short jDim=0; jDim < nDim; jDim++) {
      stress_anisotropy_tensor[iDim][jDim] = (su2double)(iDim == jDim);
    }
  }
}

void CHybrid_Isotropic_Stress::CalculateStressAnisotropy() {
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
  su2double   w, r_k;
 public:
  void SetApproxStructFunc(su2double** val_approx_struct_func);
  void SetAnisotropyWeight(su2double val_aniso_weight);

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

void CHybrid_Aniso_Q::CalculateStressAnisotropy() {
  unsigned short iDim, jDim;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      stress_anisotropy_tensor[iDim][jDim] = w*double(iDim == jDim);
      stress_anisotropy_tensor[iDim][jDim] += (1.0-w)*sqrt(3)*
          Qstar[iDim][jDim]/Qstar_norm;
    }
  }
}

inline void CHybrid_Aniso_Q::SetApproxStructFunc(su2double** val_approx_struct_func) {
  Qstar = val_approx_struct_func;
}

inline void CHybrid_Aniso_Q::SetAnisotropyWeight(su2double val_aniso_weight) {
  w = val_aniso_weight;
}

su2double** CHybrid_Aniso_Q::GetStressAnisotropyTensor() {
  return stress_anisotropy_tensor;
}







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
  virtual void SetupResolvedFlow(CGeometry* geometry,
                                 CSolver **solver_container,
                                 unsigned short iPoint) = 0;
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
  const su2double C_sf = 0.17; /*!> \brief Smagorinksy constant */
  su2double **Q_,
  **Qstar_,
  **ResolutionTensor,
  **PrimVar_Grad_i,
  **PrimVar_Grad_j;
  su2double r_k;

  /*!
   * \brief Calculates the resolution inadequacy parameter
   * @param Q - The approximate 2nd order structure function
   * @param v2 - The v2 value from Durbin's k-eps-v2-f model
   * @return The resolution inadequacy parameter
   */
  su2double CalculateRk(su2double** Q, su2double v2);

  /*!
   * \brief Calculates the approximate 2nd order structure function tensor, Qij
   * @param[in] ResolutionTensor - The resolution (metric) tensor
   * @param[in] PrimVar_Grad - Gradients of the primitive flow variables
   * @param[out] Q - The approximate 2nd order structure function tensor
   */
  void CalculateApproxStructFunc(su2double** ResolutionTensor,
                                 su2double** PrimVar_Grad,
                                 su2double** Q);

  /*!
   * \brief Calculates the anisotropy weight for the eddy viscosity tensor
   * @param[in] r_k - The resolution inadequacy parameter
   * @return The anisotropy weight for the eddy viscosity tensor
   */
  su2double CalculateAnisotropyWeight(su2double r_k);

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
   * \brief
   */
  void SetupBlendingSolver(CGeometry* geometry, CSolver **solver_container,
                           unsigned short iPoint);


  void SetupBlendingNumerics(CGeometry* geometry, CSolver **solver_container,
                             CNumerics *blending_numerics,
                             unsigned short iPoint, unsigned short jPoint);

  /**
   * \brief The turbulent stress anisotropy needs the weighting factor and the
   *        second order structure function.
   */
  void SetupStressAnisotropy(CGeometry* geometry,
                             CSolver **solver_container,
                             unsigned short iPoint);

  /**
   * \brief The resolved flow needs the blending coefficient (the energy flow
   *        parameter) and the turbulent stress anisotropy tensor.
   */
  void SetupResolvedFlow(CGeometry* geometry,
                         CSolver **solver_container,
                         CNumerics* visc_numerics,
                         unsigned short iPoint);
};


CHybrid_Mediator::CHybrid_Mediator(int nDim,
                                   CHybrid_Aniso_Q* hybrid_anisotropy)
: nDim(nDim), hybrid_anisotropy(hybrid_anisotropy) {

  //TODO: Get the Smagorinksy constant from the config file.
}

CHybrid_Mediator::~CHybrid_Mediator() {

}

void CHybrid_Mediator::SetupRANSNumerics(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         unsigned short iPoint,
                                         unsigned short jPoint) {

  // TODO: Check that this matches the muT implementation
  su2double alpha =
      solver_container[BLENDING_SOL]->node[iPoint]->GetPrimitive(0);
  rans_numerics->SetBlendingCoef(alpha);
}

void CHybrid_Mediator::SetupBlendingSolver(CGeometry* geometry,
                                           CSolver **solver_container,
                                           unsigned short iPoint) {

  /*--- Calculate and store the resolution adequacy parameter ---*/

  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double** PrimVar_Grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  su2double v2_ = solver_container[TURB_SOL]->node[iPoint]->GetPrimitive(2);

  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q_);
  su2double r_k = CalculateRk(Q_, v2_);
  solver_container[BLENDING_SOL]->node[iPoint]->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::SetupBlendingNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics *blending_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {

  /*--- Find and store turbulent length and timescales ---*/

  su2double turb_T =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbTimescale();
  su2double turb_L =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbLengthscale();

  if (turb_T <= 0) {
    cout << "Error: Turbulent timescale was " << turb_T << std::endl;
    exit(EXIT_FAILURE);
  }
  if (turb_L <= 0) {
    cout << "Error: Turbulent timescale was " << turb_L << std::endl;
    exit(EXIT_FAILURE);
  }

  blending_numerics->SetTurbTimescale(turb_T);
  blending_numerics->SetTurbLengthscale(turb_L);

  // FIXME: Do I need to pass in the resolution adequacy?
}

void CHybrid_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             unsigned short iPoint) {

}

void CHybrid_Mediator::SetupResolvedFlow(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* visc_numerics,
                                         unsigned short iPoint) {
}

su2double CHybrid_Mediator::CalculateRk(su2double** Q, su2double v2) {
  // TODO: Update this with new model.
  return 1.0;
}

void CHybrid_Mediator::CalculateApproxStructFunc(su2double** ResolutionTensor,
                                                 su2double** PrimVar_Grad,
                                                 su2double** Q) {
  unsigned int iDim, jDim, kDim, lDim, mDim;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      Q[iDim][jDim] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      for (kDim = 0; kDim < nDim; kDim++)
        for (lDim = 0; lDim < nDim; lDim++)
          for (mDim = 0; mDim < nDim; mDim++)
            Q[iDim][jDim] += ResolutionTensor[iDim][mDim]*
            PrimVar_Grad[mDim+1][kDim]*
            PrimVar_Grad[lDim+1][kDim]*
            ResolutionTensor[lDim][jDim];
}

su2double CHybrid_Mediator::CalculateAnisotropyWeight(su2double r_k) {
  return 1.0 - fmin(1.0/r_k, 1.0);
}





