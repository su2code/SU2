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
#ifdef HAVE_LAPACK
#include "mkl_lapacke.h"
#endif


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
  // TODO: Update this with new anisotropy model.
  su2double w_RANS = CalculateIsotropyWeight(resolution_adequacy);
  su2double w_LES = (1.0 - w_RANS);

  unsigned short iDim, jDim;

  // Qstar is already normalized.
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eddy_visc_anisotropy[iDim][jDim] = w_RANS*double(iDim == jDim);
      eddy_visc_anisotropy[iDim][jDim] += w_LES*sqrt(nDim)*
                                          Qstar[iDim][jDim];
    }
  }
}

inline su2double CHybrid_Aniso_Q::CalculateIsotropyWeight(su2double r_k) {
  if (r_k < 0.0) {
    cout << "ERROR: Resolution adequacy was negative! Value: " << r_k << endl;
    exit(EXIT_FAILURE);
  }
  return 1.0 - min(1.0/r_k, su2double(1.0));
}

inline void CHybrid_Aniso_Q::SetTensor(su2double** val_approx_struct_func) {
  Qstar = val_approx_struct_func;
}

inline void CHybrid_Aniso_Q::SetScalar(su2double val_r_k) {
  resolution_adequacy = val_r_k;
}



CHybrid_Mediator::CHybrid_Mediator(int nDim, CConfig* config, string filename)
 : nDim(nDim), C_sf(config->Get_Hybrid_Model_Const()) {

  /*--- Allocate the approximate structure function (used in calcs) ---*/

  Q = new su2double*[nDim];
  Qapprox = new su2double*[nDim];
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    Q[iDim] = new su2double[nDim];
    Qapprox[iDim] = new su2double[nDim];
  }

  /*--- Load the constants for mapping M to M-tilde ---*/

  if (filename == "") {
    cout << "WARNING: No file given for hybrid RANS/LES constants." << endl;
  } else {
    constants = LoadConstants(filename);
  }
}

CHybrid_Mediator::~CHybrid_Mediator() {
  for (unsigned int iDim = 0; iDim < nDim; iDim++)
    delete [] Q[iDim];
    delete [] Qapprox[iDim];
  delete [] Q;
  delete [] Qapprox;
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
  /*--- Find eigenvalues and eigenvecs for grid-based resolution tensor ---*/
  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double* ResolutionValues = geometry->node[iPoint]->GetResolutionValues();
  su2double** ResolutionVectors = geometry->node[iPoint]->GetResolutionVectors();

  /*--- Transform the approximate structure function ---*/
  su2double** PrimVar_Grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Qapprox);
  vector<vector<su2double> > zeta = BuildZeta(ResolutionValues);
  unsigned short iDim, jDim, kDim, lDim;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Q[iDim][jDim] = 0.0;
      for (kDim = 0; kDim < nDim; kDim++) {
        for (lDim = 0; lDim < nDim; lDim++) {
          Q[iDim][jDim] += zeta[iDim][kDim] * zeta[lDim][jDim] *
                           Qapprox[kDim][lDim];
        }
      }
    }
  }


  su2double min_resolved = TWO3*solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
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
  unsigned int iDim, jDim, kDim, lDim;
  /*--- Use the grid-based resolution tensor, not the anisotropy-correcting
   * resolution tensor ---*/
  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double** PrimVar_Grad =
        solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q);

  su2double total= 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      // iDim+1 to capture velocities, indexed in 1,2,3 (not 0,1,2)
      total += abs(Q[iDim][jDim]);
    }
  }
  if (total < EPS) {
    /*--- If there are negligible velocity gradients at the resolution scale,
     * set the tensor to the Kroneckor delta ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Q[iDim][jDim] = (iDim == jDim);
      }
    }
  } else {
    /*--- Compute eigenvalues and normalize by max eigenvalue ---*/
    vector<su2double> eigvalues_Q;
    vector<vector<su2double> > eigvectors_Q;
    SolveEigen(Q, eigvalues_Q, eigvectors_Q);
    vector<vector<su2double> > D(3, vector<su2double>(3));
    su2double max_eigvalue = *min_element(eigvalues_Q.begin(), eigvalues_Q.end());
    for (iDim=0; iDim < nDim; iDim++) {
      D[iDim][iDim] = eigvalues_Q[iDim]/max_eigvalue;
    }
    for (int iDim = 0; iDim < nDim; iDim++) {
      for (int jDim = 0; jDim < nDim; jDim++) {
        Q[iDim][jDim] = 0.0;
        for (int kDim = 0; kDim < nDim; kDim++) {
          for (int lDim = 0; lDim < nDim; lDim++) {
            Q[iDim][jDim] += eigvectors_Q[kDim][iDim] * D[kDim][lDim] *
                             eigvectors_Q[lDim][jDim];
          }
        }
      }
    }
  }

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

  su2double* alpha_i = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  su2double* alpha_j = solver_container[HYBRID_SOL]->node[jPoint]->GetSolution();
  visc_numerics->SetHybridParameter(alpha_i, alpha_j);

  /*--- Pass the stress anisotropy tensor to the resolved flow ---*/

  su2double** aniso_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscAnisotropy();
  su2double** aniso_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscAnisotropy();
  visc_numerics->SetEddyViscAnisotropy(aniso_i, aniso_j);
}

su2double CHybrid_Mediator::CalculateRk(su2double** Q,
                                        su2double min_resolved) {
  unsigned int iDim, jDim;
  su2double total = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      total += abs(Q[iDim][jDim]);
    }
  }
  if (total < EPS) {
    /*--- If the velocity gradients are negligible, return 0 ---*/
    return 0.0;
  } else {
    /*--- Compare resolved and unresolved turbulence ---*/
    vector<su2double> eigvalues_Q;
    vector<vector<su2double> > eigvectors_Q;
    SolveEigen(Q, eigvalues_Q, eigvectors_Q);
    su2double max_eigenvalue_Q = *max_element(eigvalues_Q.begin(),
                                              eigvalues_Q.end());
    su2double max_unresolved = C_sf*TWO3*max_eigenvalue_Q;
    return max_unresolved / min_resolved;
  }
}

vector<vector<su2double> > CHybrid_Mediator::LoadConstants(string filename) {
  vector<vector<su2double> > output;
  output.resize(nDim);
  ifstream file;
  for (int i=0; i<nDim; i++) {
    stringstream ss;
    ss << filename << i << ".dat";
    string fullname = ss.str();
    file.open(fullname.c_str());
    if (file.is_open()) {
      su2double value;
      while (file >> value) {
        output[i].push_back(value);
      }
      file.close();
    } else {
      cout << "ERROR: Could not open the hybrid constants file." << endl;
      cout << "       Tried reading file " << fullname << endl;
      exit(EXIT_FAILURE);
    }
  }
  return output;
};

vector<su2double> CHybrid_Mediator::GetEigValues_Q(vector<su2double> eigvalues_M) {
  su2double dnorm = *min_element(eigvalues_M.begin(), eigvalues_M.end());

  /*--- Normalize eigenvalues ---*/
  su2double a, b;
  if (eigvalues_M[0] == dnorm) {
    a = eigvalues_M[1]/dnorm;
    b = eigvalues_M[2]/dnorm;
  } else if (eigvalues_M[1] == dnorm) {
    a = eigvalues_M[0]/dnorm;
    b = eigvalues_M[2]/dnorm;
  } else {
    a = eigvalues_M[0]/dnorm;
    b = eigvalues_M[1]/dnorm;
  }

  /*--- Convert to cylindrical coordinates ---*/
  su2double r = sqrt(a*a + b*b);
  su2double theta = acos(max(a,b)/r);

  /*--- Convert to more convenient log coordinates ---*/
  su2double x = log(sin(2*theta));
  su2double y = log(r);

  vector<su2double> g(3);
  for (int iDim = 0; iDim < nDim; iDim++) {
    g[iDim] = constants[iDim][0];
    g[iDim] += constants[iDim][1]*x;
    g[iDim] += constants[iDim][2]*y;
    g[iDim] += constants[iDim][3]*x*x;
    g[iDim] += constants[iDim][4]*x*y;
    g[iDim] += constants[iDim][5]*y*y;
    g[iDim] += constants[iDim][6]*x*x*x;
    g[iDim] += constants[iDim][7]*x*x*y;
    g[iDim] += constants[iDim][8]*x*y*y;
    g[iDim] += constants[iDim][9]*y*y*y;
    g[iDim] += constants[iDim][10]*x*x*x*x;
    g[iDim] += constants[iDim][11]*x*x*x*y;
    g[iDim] += constants[iDim][12]*x*x*y*y;
    g[iDim] += constants[iDim][13]*x*y*y*y;
    g[iDim] += constants[iDim][14]*y*y*y*y;
  }

  vector<su2double> eigvalues_Q(3);
  for (int iDim = 0; iDim < nDim; iDim++) {
    eigvalues_Q[iDim] = g[0] + g[1]*log(eigvalues_M[iDim]/dnorm) +
                        g[2]*pow(log(eigvalues_M[iDim]/dnorm),2);
    eigvalues_Q[iDim] = exp(eigvalues_Q[iDim]);
  }

  return eigvalues_Q;
}

vector<vector<su2double> > CHybrid_Mediator::BuildZeta(su2double* values_M) {

#ifdef HAVE_LAPACK
  unsigned short iDim;

  vector<su2double> eigvalues_M(nDim);
  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues_M[iDim] = values_M[iDim];

  /*--- Solve for the modified resolution tensor  ---*/
  vector<su2double> eigvalues_Mtilde = GetEigValues_Zeta(eigvalues_M);
  vector<vector<su2double> > D(3, vector<su2double>(3));
  for (iDim = 0; iDim < nDim; iDim++) {
    zeta[iDim][iDim] = eigvalues_Zeta[iDim];
  }
  return zeta;
#else
  cout << "Eigensolver without LAPACK not implemented; please use LAPACK." << endl;
  exit(EXIT_FAILURE);
#endif

}

vector<su2double> CHybrid_Mediator::GetEigValues_Zeta(vector<su2double> eigvalues_M) {
  su2double C_zeta = 0.46120398645; // Overall scaling constant

  /*--- Find the minimum eigenvalue ---*/

  su2double dnorm = *min_element(eigvalues_M.begin(), eigvalues_M.end());

  /*--- Get eigenvalues to the normalized gradient-gradien tensor ---*/

  vector<su2double> eigvalues_Q = GetEigValues_Q(eigvalues_M);

  /*--- Use eigenvalues from M and G to find eigenvalues for modified M ---*/

  vector<su2double> eigvalues_zeta(3);
  for (int iDim = 0; iDim < nDim; iDim++) {
    eigvalues_zeta[iDim] = C_zeta*pow((eigvalues_M[iDim]/dnorm),1.0/3)*
                           pow(eigvalues_Q[iDim],-0.5);
  }

  return eigvalues_zeta;
}

void CHybrid_Mediator::SolveEigen(su2double** M,
                                  vector<su2double> &eigvalues,
                                  vector<vector<su2double> > &eigvectors) {
unsigned short iDim, jDim;

#ifndef NDEBUG
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (M[iDim][jDim] != M[iDim][jDim]) {
        cout << "ERROR: SolveEigen received a matrix with NaN!" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  su2double sum = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      sum += abs(M[iDim][jDim]);
    }
  }
  if (sum < EPS) {
    cout << "ERROR: SolveEigen received an empty matrix!" << endl;
    exit(EXIT_FAILURE);
  }
#endif


  eigvalues.resize(nDim);
  eigvectors.resize(nDim, std::vector<su2double>(nDim));

#ifdef HAVE_LAPACK
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim< nDim; jDim++) {
      mat[iDim*nDim+jDim] = M[iDim][jDim];
    }
  }

  /*--- Call LAPACK routines ---*/
  info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', nDim, mat, nDim, wr, wi, vl, nDim, vr, nDim);
  if (info != 0) {
    cout << "ERROR: The solver failed to compute eigenvalues." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Check the values ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    if (wr[iDim] < 0.0) {
      if (wr[iDim] > -1e-6) {
        wr[iDim] = 0.0;
      } else {
        cout << "ERROR: The solver return a large negative eigenvalue!" << endl;
        cout << "    Eigenvalues:" << endl;
        cout << "    [" << endl;
        cout << wr[0] << ", " << endl;
        cout << wr[1] << ", " << endl;
        cout << wr[2] << "]" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  /*--- Rewrite arrays to eigenvalues/vectors output ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues[iDim] = wr[iDim];

  for (iDim = 0; iDim < nDim; iDim++) {
    su2double norm = 0.0;
    for (jDim = 0; jDim < nDim; jDim++)
      norm += vr[iDim*nDim+jDim]*vr[iDim*nDim+jDim];
    norm = sqrt(norm);
    for (jDim = 0; jDim < nDim; jDim++) {
      eigvectors[iDim][jDim] = vr[iDim*nDim+jDim]/norm;
    }
  }
#else
  cout << "Eigensolver without LAPACK not implemented; please use LAPACK." << endl;
  exit(EXIT_FAILURE);
#endif

}

vector<vector<su2double> > CHybrid_Mediator::GetConstants() {
  return constants;
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
