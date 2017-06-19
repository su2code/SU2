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

CHybrid_Visc_Anisotropy::CHybrid_Visc_Anisotropy(unsigned short nDim)
    : nDim(nDim) {
  eddy_visc_anisotropy = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    eddy_visc_anisotropy[iDim] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++)
      eddy_visc_anisotropy[iDim][jDim] = 0.0;
  }
}

CHybrid_Visc_Anisotropy::~CHybrid_Visc_Anisotropy() {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] eddy_visc_anisotropy[iDim];
  delete [] eddy_visc_anisotropy;
}

su2double** CHybrid_Visc_Anisotropy::GetViscAnisotropy() {
  return eddy_visc_anisotropy;
}

CHybrid_Isotropic_Visc::CHybrid_Isotropic_Visc(unsigned short nDim)
: CHybrid_Visc_Anisotropy(nDim) {
  for (unsigned short iDim=0; iDim < nDim; iDim++) {
    for (unsigned short jDim=0; jDim < nDim; jDim++) {
      eddy_visc_anisotropy[iDim][jDim] = (su2double)(iDim == jDim);
    }
  }
}

void CHybrid_Isotropic_Visc::CalculateViscAnisotropy() {
};

CHybrid_Aniso_Q::CHybrid_Aniso_Q(unsigned short nDim)
  : CHybrid_Visc_Anisotropy(nDim) {
}

void CHybrid_Aniso_Q::CalculateViscAnisotropy() {
  su2double w_RANS = CalculateIsotropyWeight(resolution_adequacy);

  Qstar_norm = (Qstar[0][0] + Qstar[1][1] + Qstar[2][2])/3.0;

  unsigned short iDim, jDim;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eddy_visc_anisotropy[iDim][jDim] = w_RANS*double(iDim == jDim);
      eddy_visc_anisotropy[iDim][jDim] += (1.0-w_RANS)*sqrt(3)*
          Qstar[iDim][jDim]/Qstar_norm;
    }
  }
}

inline su2double CHybrid_Aniso_Q::CalculateIsotropyWeight(su2double r_k) {
  return 1.0 - fmin(double(1.0/r_k), 1.0);
}

inline void CHybrid_Aniso_Q::SetTensor(su2double** val_approx_struct_func) {
  Qstar = val_approx_struct_func;
}

inline void CHybrid_Aniso_Q::SetScalar(su2double val_r_k) {
  resolution_adequacy = val_r_k;
}



CHybrid_Mediator::CHybrid_Mediator(int nDim, CConfig* config)
 : nDim(nDim), C_sf(config->Get_Hybrid_Model_Const()) {

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
      solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  /*--- Since this is a source term, we don't need a second point ---*/
  rans_numerics->SetHybridParameter(alpha, NULL);
}

void CHybrid_Mediator::SetupHybridParamSolver(CGeometry* geometry,
                                              CSolver **solver_container,
                                              unsigned short iPoint) {

  /*--- Calculate and store the resolution adequacy parameter ---*/

  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double** PrimVar_Grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  // su2double min_resolved  = solver_container[TURB_SOL]->node[iPoint]->GetPrimitive(2);
  su2double min_resolved = TWO3*solver_container[TURB_SOL]->node[iPoint]->GetPrimitive(2);

  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q);
  su2double r_k = CalculateRk(Q, min_resolved);
  solver_container[HYBRID_SOL]->node[iPoint]->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::SetupHybridParamNumerics(CGeometry* geometry,
                                                CSolver **solver_container,
                                                CNumerics *hybrid_param_numerics,
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

  hybrid_param_numerics->SetTurbTimescale(turb_T);
  hybrid_param_numerics->SetTurbLengthscale(turb_L);

  /*--- Pass resolution adequacy into the numerics object ---*/

  su2double r_k = solver_container[HYBRID_SOL]->node[iPoint]->GetResolutionAdequacy();
  hybrid_param_numerics->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                                             unsigned short iPoint) {

  /*--- Find Approximate Structure Function ---*/

  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double** PrimVar_Grad =
        solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q);
  hybrid_anisotropy->SetTensor(Q);

  /*--- Retrieve and pass along the resolution adequacy parameter ---*/

  su2double r_k = solver_container[HYBRID_SOL]->node[iPoint]->GetResolutionAdequacy();
  hybrid_anisotropy->SetScalar(r_k);
}

void CHybrid_Mediator::SetupResolvedFlowNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {

  /*--- Pass alpha to the resolved flow ---*/

  su2double* alpha_i = solver_container[HYBRID_SOL]->node[iPoint]->GetPrimitive();
  su2double* alpha_j = solver_container[HYBRID_SOL]->node[jPoint]->GetPrimitive();
  visc_numerics->SetHybridParameter(alpha_i, alpha_j);

  /*--- Pass the stress anisotropy tensor to the resolved flow ---*/

  su2double** aniso_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscAnisotropy();
  su2double** aniso_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscAnisotropy();
  visc_numerics->SetEddyViscAnisotropy(aniso_i, aniso_j);
}

su2double CHybrid_Mediator::CalculateRk(su2double** Q,
                                        su2double min_resolved) {
  // TODO: Implement function here to find max eigienvalue of Q.
  su2double max_eigenvalue_Q = 1.0;
  su2double max_unresolved = C_sf*TWO3*max_eigenvalue_Q;
  return max_unresolved / min_resolved;
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



CHybrid_Dummy_Mediator::CHybrid_Dummy_Mediator(int nDim, CConfig* config) : nDim(nDim) {

  /*--- Set the default value of the hybrid parameter ---*/
  dummy_alpha = new su2double[1];
  dummy_alpha[0] = 1.0;
}

CHybrid_Dummy_Mediator::~CHybrid_Dummy_Mediator() {
  delete [] dummy_alpha;
}

void CHybrid_Dummy_Mediator::SetupRANSNumerics(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         unsigned short iPoint,
                                         unsigned short jPoint) {
  // This next line is just here for testing purposes.
  su2double* alpha =
      solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  rans_numerics->SetHybridParameter(dummy_alpha, dummy_alpha);
}

void CHybrid_Dummy_Mediator::SetupHybridParamSolver(CGeometry* geometry,
                                           CSolver **solver_container,
                                           unsigned short iPoint) {
}

void CHybrid_Dummy_Mediator::SetupHybridParamNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics *hybrid_param_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {
}

void CHybrid_Dummy_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                                             unsigned short iPoint) {
}

void CHybrid_Dummy_Mediator::SetupResolvedFlowNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {

  /*--- Pass alpha to the resolved flow ---*/

  // This next two lines are just here for testing purposes.
  su2double* alpha_i = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  su2double* alpha_j = solver_container[HYBRID_SOL]->node[jPoint]->GetSolution();
  visc_numerics->SetHybridParameter(dummy_alpha, dummy_alpha);

  /*--- Pass the stress anisotropy tensor to the resolved flow ---*/

  su2double** aniso_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscAnisotropy();
  su2double** aniso_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscAnisotropy();
  visc_numerics->SetEddyViscAnisotropy(aniso_i, aniso_j);
}
