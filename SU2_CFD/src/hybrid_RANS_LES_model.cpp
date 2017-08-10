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
#include "mkl.h"
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
}

CHybrid_Aniso_Q::CHybrid_Aniso_Q(unsigned short nDim)
  : rans_weight(1.0), CHybrid_Visc_Anisotropy(nDim) {
}

void CHybrid_Aniso_Q::CalculateViscAnisotropy() {
#ifndef NDEBUG
  if (rans_weight < 0 || rans_weight > 1) {
    cout << "ERROR: Isotropic weighting in hybrid RANS/LES was not in range [0,1]" << endl;
    cout << "       weight = " << rans_weight << endl;
    exit(EXIT_FAILURE);
  }
#endif
  su2double LES_weight = (1.0 - rans_weight);

  unsigned short iDim, jDim;

  // Qstar is already normalized.
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eddy_visc_anisotropy[iDim][jDim] = rans_weight*double(iDim == jDim);
      eddy_visc_anisotropy[iDim][jDim] += LES_weight*sqrt(nDim) *
                                          Qstar[iDim][jDim];
    }
  }
}


inline void CHybrid_Aniso_Q::SetTensor(su2double** val_approx_struct_func) {
  Qstar = val_approx_struct_func;
}

inline void CHybrid_Aniso_Q::SetScalars(vector<su2double> val_scalars) {
#ifndef NDEBUG
  if (val_scalars.size() != 1) {
    cout << "ERROR: Improper number of scalars passed to anisotropy model!" << endl;
    cout << "       Expected: 1" << endl;
    cout << "       Found:    " << val_scalars.size() << endl;
    exit(EXIT_FAILURE);
  }
#endif
  rans_weight = val_scalars[0];
}



CHybrid_Mediator::CHybrid_Mediator(int nDim, CConfig* config, string filename)
 : nDim(nDim), C_sf(config->Get_Hybrid_Model_Const()) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Allocate the approximate structure function (used in calcs) ---*/

  Q = new su2double*[nDim];
  Qapprox = new su2double*[nDim];
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    Q[iDim] = new su2double[nDim];
    Qapprox[iDim] = new su2double[nDim];
  }

  /*--- Load the constants for mapping M to M-tilde ---*/

  if (filename == "") {
    if (rank == MASTER_NODE) {
      cout << "WARNING: No file given for hybrid RANS/LES constants." << endl;
      cout << "         Default (hardcoded) values used." << endl;
    }
    constants.resize(3, vector<su2double>(15));
    constants[0][0]  = -1.3441989500430;
    constants[0][1]  = -1.2997962353961;
    constants[0][2]  =  0.3624921526914;
    constants[0][3]  = -0.0114257653000;
    constants[0][4]  = -0.0002359408737;
    constants[0][5]  =  0.2314555700009;
    constants[0][6]  =  0.0016951253508;
    constants[0][7]  = -0.0002866338640;
    constants[0][8]  = -0.0009095648353;
    constants[0][9]  =  0.0617742005401;
    constants[0][10] = -0.0000918936400;
    constants[0][11] =  0.0000544349543;
    constants[0][12] =  0.0001553623566;
    constants[0][13] =  0.0001159449506;
    constants[0][14] =  0.0059904950363;

    constants[1][0]  =  1.9527214390792;
    constants[1][1]  =  0.0387097819862;
    constants[1][2]  = -0.0035236439484;
    constants[1][3]  = -0.0090288063753;
    constants[1][4]  =  0.0074191857589;
    constants[1][5]  =  0.0012129089188;
    constants[1][6]  =  0.0010376720357;
    constants[1][7]  = -0.0011418423755;
    constants[1][8]  =  0.0008778876148;
    constants[1][9]  =  0.0006097328257;
    constants[1][10] = -0.0000471680743;
    constants[1][11] =  0.0000668786855;
    constants[1][12] = -0.0000477449682;
    constants[1][13] =  0.0000549683602;
    constants[1][14] =  0.0000759797725;

    constants[2][0]  = -0.0149586909340;
    constants[2][1]  =  0.0021107470732;
    constants[2][2]  = -0.0057740374642;
    constants[2][3]  =  0.0002190519772;
    constants[2][4]  =  0.0006652247252;
    constants[2][5]  = -0.0017709687027;
    constants[2][6]  = -0.0000860428032;
    constants[2][7]  =  0.0000032821371;
    constants[2][8]  =  0.0001955198753;
    constants[2][9]  = -0.0002762253983;
    constants[2][10] =  0.0000061997774;
    constants[2][11] = -0.0000041785763;
    constants[2][12] = -0.0000092567404;
    constants[2][13] =  0.0000125295721;
    constants[2][14] = -0.0000180552959;
  
  } else {
    constants = LoadConstants(filename);
  }

  /*--- Calculate scaling constants so that zeta -> kroneckor delta for
   * isotropic cells ---*/
  vector<su2double> temp_values = GetEigValues_Q(vector<su2double>(3, 1.0));
  C_zeta = pow(temp_values[0],0.5);
}

CHybrid_Mediator::~CHybrid_Mediator() {
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    delete [] Q[iDim];
    delete [] Qapprox[iDim];
  }
  delete [] Q;
  delete [] Qapprox;

#ifdef HAVE_LAPACK
  mkl_free_buffers();
#endif
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
  unsigned short iDim, jDim, kDim, lDim;


  /*--- Find eigenvalues and eigenvecs for grid-based resolution tensor ---*/
  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double* ResolutionValues = geometry->node[iPoint]->GetResolutionValues();
  su2double** ResolutionVectors = geometry->node[iPoint]->GetResolutionVectors();

  /*--- Transform the approximate structure function ---*/
  su2double** PrimVar_Grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Qapprox);
  vector<vector<su2double> > zeta = BuildZeta(ResolutionValues);
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Q[iDim][jDim] = 0.0;
      for (kDim = 0; kDim < nDim; kDim++) {
        for (lDim = 0; lDim < nDim; lDim++) {
          // Now zeta*Q*zeta, not Q
          Q[iDim][jDim] += zeta[iDim][kDim] * zeta[lDim][jDim] *
                           Qapprox[kDim][lDim];
        }
      }
    }
  }

  /*--- Find eigenvalues and eigenvectors ---*/
  su2double r_k, w_rans;
  su2double total_vel_differences = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      total_vel_differences += abs(Q[iDim][jDim]);
    }
  }
  if (geometry->node[iPoint]->GetSolidBoundary()) {
    /*--- At solid walls, force the RANS weight to 1 and the
     * resolution adequacy to 0 (no resolved fluctuations) ---*/
    r_k = 0.0;
    w_rans = 1.0;

  } else if (total_vel_differences > EPS) {
    /*--- Only calculate r_k, w_rans if there are resolved velocity differences
     * at resolution scale.  Otherwise, eigenvector calculation is arbitrary */

    /*--- Calculate eigenvectors and eigenvalues of zeta*Q*zeta ---*/
    vector<su2double> eigvalues_zQz;
    vector<vector<su2double> > eigvectors_zQz;
    SolveEigen(Q, eigvalues_zQz, eigvectors_zQz);
    std::vector<su2double>::iterator iter;
    iter = max_element(eigvalues_zQz.begin(), eigvalues_zQz.end());
    unsigned short max_index = distance(eigvalues_zQz.begin(), iter);
    vector<su2double> max_eigenvalue_direction = eigvectors_zQz[max_index];

    /*---Find the largest product of resolved fluctuations at the cutoff---*/
    // TODO: Incorporate anisotropy ratio here (max-to-min)
    su2double aniso_ratio = 1.0;
    su2double C_kQ = 8.0;
    su2double max_resolved = aniso_ratio*C_kQ*C_sf*TWO3*eigvalues_zQz[max_index];

    /*--- Find the smallest product of unresolved fluctuations at the cutoff ---*/
    // TODO: Make this more generalized
    su2double C_mu = 0.22;
    su2double TurbT = solver_container[TURB_SOL]->node[iPoint]->GetTurbTimescale();
    su2double k = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    su2double omega = solver_container[TURB_SOL]->node[iPoint]->GetSolution(1);
    su2double min_unresolved = TurbT*k*omega/C_mu;

    /*--- Calculate the resolution adequacy parameter ---*/
    r_k = max_resolved/min_unresolved;

    /*--- Find the dissipation ratio ---*/
    su2double C_eps = 0.03125;
    su2double TurbL = solver_container[TURB_SOL]->node[iPoint]->GetTurbLengthscale();
    su2double d_max = GetProjResolution(ResolutionTensor,
                                        eigvectors_zQz[max_index]);
    su2double r_eps = C_eps * pow(r_k, 1.5) * TurbL / d_max;

    /*--- Calculate the RANS weight ---*/
    w_rans = max(tanh(r_eps - 1), 0.0);

  } else {
    /*--- If the velocity differences at resolution scale are negligible,
     * set max_eigenvalue = 0 ---*/
    r_k = 0.0;
    w_rans = 0.0;
  }

  solver_container[HYBRID_SOL]->node[iPoint]->SetResolutionAdequacy(r_k);
  solver_container[HYBRID_SOL]->node[iPoint]->SetRANSWeight(w_rans);

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

#ifndef NDEBUG
  if (turb_T <= 0) {
    cout << "Error: Turbulent timescale was " << turb_T << std::endl;
    exit(EXIT_FAILURE);
  }
  if (turb_L <= 0) {
    cout << "Error: Turbulent timescale was " << turb_L << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  hybrid_param_numerics->SetTurbTimescale(turb_T);
  hybrid_param_numerics->SetTurbLengthscale(turb_L);

  /*--- Pass resolution adequacy into the numerics object ---*/

  su2double r_k = solver_container[HYBRID_SOL]->node[iPoint]->GetResolutionAdequacy();
  hybrid_param_numerics->SetResolutionAdequacy(r_k);

  /*--- Pass RANS weight into the numerics object ---*/

  su2double w_rans = solver_container[HYBRID_SOL]->node[iPoint]->GetRANSWeight();
  hybrid_param_numerics->SetRANSWeight(w_rans);
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
      total += abs(Q[iDim][jDim]);
    }
  }
  if (total < 9*EPS) {
    /*--- If there are negligible velocity differences at the resolution scale,
     * set the tensor to the Kroneckor delta ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Q[iDim][jDim] = (iDim == jDim);
      }
    }
  } else {
    /*--- Normalize the approximate structure function tensor ---*/
    su2double norm = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        norm += Q[iDim][jDim]*Q[jDim][iDim];
      }
    }
    norm = sqrt(norm);
    for (iDim=0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Q[iDim][jDim] /= norm;
      }
    }
  }

#ifndef NDEBUG
  su2double trace = Q[0][0] + Q[1][1] + Q[2][2];
  if (trace > 1e7 || trace < 1e-7) {
      cout << "ERROR: Trace of the stress anisotropy was unusual." << endl;
      cout << "       Trace: " << trace << endl;
      exit(EXIT_FAILURE);
  }
#endif

  hybrid_anisotropy->SetTensor(Q);

  /*--- Retrieve and pass along the resolution adequacy parameter ---*/

  vector<su2double> scalars;
  scalars.push_back(solver_container[HYBRID_SOL]->node[iPoint]->GetRANSWeight());
  hybrid_anisotropy->SetScalars(scalars);
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

su2double CHybrid_Mediator::GetProjResolution(su2double** resolution_tensor,
                                              vector<su2double> direction) {

#ifndef NDEBUG
  su2double magnitude_squared = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    magnitude_squared += direction[iDim]*direction[iDim];
  if (abs(magnitude_squared - 1.0) > 1e-7) {
    cout << "ERROR: The unit vector int the projected resolution calc had a ";
    cout << "mangitude greater than 1!" << endl;
    cout << "    Magnitude: " << sqrt(magnitude_squared) << endl;
    exit(EXIT_FAILURE);
  }
#endif

  su2double temp, result = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    temp = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      temp += resolution_tensor[iDim][jDim]*direction[jDim];
    }
    result += temp*temp;
  }
  return sqrt(result);
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

  vector<su2double> eigvalues_M;
  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues_M.push_back(values_M[iDim]);

  /*--- Solve for the modified resolution tensor  ---*/
  vector<su2double> eigvalues_Zeta = GetEigValues_Zeta(eigvalues_M);
  vector<vector<su2double> > zeta(3, vector<su2double>(3));
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
  /*--- Find the minimum eigenvalue ---*/

  su2double dnorm = *min_element(eigvalues_M.begin(), eigvalues_M.end());

  /*--- Get eigenvalues to the normalized gradient-gradien tensor ---*/

  vector<su2double> eigvalues_Q = GetEigValues_Q(eigvalues_M);
#ifndef NDEBUG
  for (int iDim = 0; iDim < nDim; iDim++) {
    if (eigvalues_Q[iDim] != eigvalues_Q[iDim]) {
      cout << "ERROR: At least one computed eigenvalue was NaN!" << endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

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

  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (iDim == jDim) continue;
      if (abs(M[iDim][jDim] - M[jDim][iDim]) > 1e-10 &&
          abs(M[iDim][jDim] - M[jDim][iDim])/abs(M[iDim][jDim]) > 1e-6) {
        cout << "ERROR: SolveEigen received a non-symmetric matrix!" << endl;
        cout << "    The difference between elements" << endl;
        cout << "      [" << iDim << ", " << jDim << "] and [" << jDim << ", " << iDim << "]" << endl;
        cout << "      was: " << M[iDim][jDim] - M[jDim][iDim] << endl;
        cout << "    Matrix:" << endl;
        cout << "      [[" << M[0][0] << ", " << M[0][1] << ", " << M[0][2] << "]" << endl;
        cout << "       [" << M[1][0] << ", " << M[1][1] << ", " << M[1][2] << "]" << endl;
        cout << "       [" << M[2][0] << ", " << M[2][1] << ", " << M[2][2] << "]]" << endl;
        exit(EXIT_FAILURE);
      }
    }
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

  info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nDim, mat, nDim, eigval);
  if (info != 0) {
    cout << "ERROR: The solver failed to compute eigenvalues." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Rewrite arrays to eigenvalues output ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues[iDim] = eigval[iDim];

  /*--- Check the values ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    if (eigval[iDim] < 0.0 && eigval[iDim] > -1e-6) {
        eigvalues[iDim] = 0.0;
    } else if (eigval[iDim] < -1e-6) {
      cout << "ERROR: The solver returned a large negative eigenvalue!" << endl;
      cout << "    Eigenvalues: [";
      cout << eigval[0] << ", ";
      cout << eigval[1] << ", ";
      cout << eigval[2] << "]" << endl;
      exit(EXIT_FAILURE);
    }
  }

  su2double max_val = max(eigval[0], max(eigval[1], eigval[2]));
  for (iDim = 0; iDim < nDim; iDim++) {
    if (eigvalues[iDim] > 0 && log10(max_val/eigvalues[iDim]) > 10) {
      // If the condition number is bad, give up
      eigvalues[iDim] = 0.0;
    }
  }

  /*--- Normalize the eigenvectors by the L2 norm of each vector ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    su2double norm = 0.0;
    for (jDim = 0; jDim < nDim; jDim++)
      norm += mat[iDim*nDim+jDim]*mat[iDim*nDim+jDim];
    norm = sqrt(norm);
    for (jDim = 0; jDim < nDim; jDim++) {
      eigvectors[iDim][jDim] = mat[iDim*nDim+jDim]/norm;
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
