/*!
 * \file hybrid_RANS_LES_model.cpp
 * \brief Describes the hybrid RANS/LES models
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

#include "../include/hybrid_RANS_LES_model.hpp"
#include "../include/numerics_structure.hpp"
#include "../include/solver_structure.hpp"

CHybrid_SGS_Anisotropy::CHybrid_SGS_Anisotropy(unsigned short nDim)
    : nDim(nDim) {
  stress_anisotropy_tensor = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    stress_anisotropy_tensor[iDim] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++)
      stress_anisotropy_tensor[iDim][jDim] = 0.0;
  }
}

CHybrid_SGS_Anisotropy::~CHybrid_SGS_Anisotropy() {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] stress_anisotropy_tensor[iDim];
  delete [] stress_anisotropy_tensor;
}

su2double** CHybrid_SGS_Anisotropy::GetStressAnisotropy() {
  return stress_anisotropy_tensor;
}

CHybrid_Isotropic_Stress::CHybrid_Isotropic_Stress(unsigned short nDim)
: CHybrid_SGS_Anisotropy(nDim) {
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

CHybrid_Aniso_Q::CHybrid_Aniso_Q(unsigned short nDim)
  : CHybrid_SGS_Anisotropy(nDim) {
}

void CHybrid_Aniso_Q::CalculateStressAnisotropy() {
  su2double w_RANS = CalculateIsotropyWeight(resolution_adequacy);

  // FIXME: How to get Qstar_norm?
  Qstar_norm = Qstar[0][0] + Qstar[1][1] + Qstar[2][2];

  unsigned short iDim, jDim;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      stress_anisotropy_tensor[iDim][jDim] = w_RANS*double(iDim == jDim);
      stress_anisotropy_tensor[iDim][jDim] += (1.0-w_RANS)*sqrt(3)*
          Qstar[iDim][jDim]/Qstar_norm;
    }
  }
}

inline su2double CHybrid_Aniso_Q::CalculateIsotropyWeight(su2double r_k) {
  return 1.0 - fmin(1.0/r_k, 1.0);
}

inline void CHybrid_Aniso_Q::SetApproxStructFunc(su2double** val_approx_struct_func) {
  Qstar = val_approx_struct_func;
}

inline void CHybrid_Aniso_Q::SetResolutionAdequacy(su2double val_r_k) {
  resolution_adequacy = val_r_k;
}

CHybrid_Mediator::CHybrid_Mediator(int nDim,
                                   CHybrid_Aniso_Q* hybrid_anisotropy)
: nDim(nDim), hybrid_anisotropy(hybrid_anisotropy), C_sf(0.17) {
  
  //TODO: Get the Smagorinksy constant from the config file.

  /*--- Allocate the approximate structure function (used in calcs) ---*/
  Q = new su2double*[nDim];
  for (unsigned int iDim = 0; iDim < nDim; iDim++)
    Q[iDim] = new su2double[nDim];
}

CHybrid_Mediator::~CHybrid_Mediator() {
  for (unsigned int iDim = 0; iDim < nDim; iDim++)
    delete [] Q[iDim];
  delete [] Q;
}

void CHybrid_Mediator::SetupRANSNumerics(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         unsigned short iPoint,
                                         unsigned short jPoint) {
  su2double* alpha =
      solver_container[BLEND_SOL]->node[iPoint]->GetPrimitive();
  // TODO: Check what other source term functions do for Set/Get
  rans_numerics->SetBlendingCoef(alpha, alpha);
}

void CHybrid_Mediator::SetupBlendingSolver(CGeometry* geometry,
                                           CSolver **solver_container,
                                           unsigned short iPoint) {

  /*--- Calculate and store the resolution adequacy parameter ---*/

  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double** PrimVar_Grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  su2double v2_ = solver_container[TURB_SOL]->node[iPoint]->GetPrimitive(2);

  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q);
  su2double r_k = CalculateRk(Q, v2_);
  solver_container[BLEND_SOL]->node[iPoint]->SetResolutionAdequacy(r_k);
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

  /*--- Pass resolution adequacy into the numerics object ---*/

  su2double r_k = solver_container[BLEND_SOL]->node[iPoint]->GetSolution(0);
  blending_numerics->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CHybrid_SGS_Anisotropy* hybrid_anisotropy,
                                             unsigned short iPoint) {

  /*--- Find Approximate Structure Function ---*/

  // FIXME: Difference between Q* and Q?
  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double** PrimVar_Grad =
        solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q);
  hybrid_anisotropy->SetApproxStructFunc(Q);

  /*--- Retrieve and pass along the resolution adequacy parameter ---*/

  su2double r_k = solver_container[BLEND_SOL]->node[iPoint]->GetResolutionAdequacy();
  hybrid_anisotropy->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::SetupMeanFlowNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {

  /*--- Pass alpha to the mean flow ---*/

  su2double* alpha_i = solver_container[BLEND_SOL]->node[iPoint]->GetSolution();
  su2double* alpha_j = solver_container[BLEND_SOL]->node[jPoint]->GetSolution();
  visc_numerics->SetBlendingCoef(alpha_i, alpha_j);

  /*--- Pass the stress anisotropy tensor to the mean flow ---*/

  su2double** aniso_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscAnisotropy();
  su2double** aniso_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscAnisotropy();
  visc_numerics->SetEddyViscAnisotropy(aniso_i, aniso_j);
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
