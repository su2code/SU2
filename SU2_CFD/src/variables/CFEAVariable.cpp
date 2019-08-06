/*!
 * \file CFEAVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
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

#include "../../include/variables/CFEAVariable.hpp"


//CFEAVariable::CFEAVariable(su2double *val_fea, Idx_t val_nDim, Idx_t val_nvar,
//                           CConfig *config) : CVariable(val_nDim, val_nvar, config) {
//
//  Idx_t iVar;
//  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
//  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
//  bool incremental_load = config->GetIncrementalLoad();
//  bool prestretch_fem = config->GetPrestretch();    // Structure is prestretched
//
//  bool discrete_adjoint = config->GetDiscrete_Adjoint();
//
//  bool refgeom = config->GetRefGeom();        // Reference geometry needs to be stored
//
//  bool dynamic_analysis = (config->GetDynamic_Analysis() == DYNAMIC);
//  bool fsi_analysis = config->GetFSI_Simulation();
//
//  VonMises_Stress       = 0.0;
//
//  Stress                = NULL;    // Nodal stress (for output purposes)
//  Residual_Ext_Body     = NULL;    // Residual component due to body forces
//
//  Solution_Vel          = NULL;    // Velocity at the node at time t+dt
//  Solution_Vel_time_n   = NULL;    // Velocity at the node at time t
//
//  Solution_Accel        = NULL;    // Acceleration at the node at time t+dt
//  Solution_Accel_time_n = NULL;    // Acceleration at the node at time t
//
//  Solution_Pred         = NULL;    // Predictor of the solution at the current subiteration
//  Solution_Pred_Old     = NULL;    // Predictor of the solution at the previous subiteration
//
//  Prestretch            = NULL;    // Prestretch geometry
//  Reference_Geometry    = NULL;    // Reference geometry for optimization purposes
//
//  Solution_BGS_k        = NULL;    // Old solution stored to check convergence in the BGS loop
//
//  if      (nDim == 2) Stress = new su2double [3];
//  else if (nDim == 3) Stress = new su2double [6];
//
//  /*--- Initialization of variables ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    Solution[iVar] = val_fea[iVar];
//  }
//
//  if (dynamic_analysis) {
//    Solution_Vel          =  new su2double [nVar];
//    Solution_Vel_time_n   =  new su2double [nVar];
//    Solution_Accel        =  new su2double [nVar];
//    Solution_Accel_time_n =  new su2double [nVar];
//    for (iVar = 0; iVar < nVar; iVar++) {
//      Solution_Vel[iVar]          = val_fea[iVar+nVar];
//      Solution_Vel_time_n[iVar]   = val_fea[iVar+nVar];
//      Solution_Accel[iVar]        = val_fea[iVar+2*nVar];
//      Solution_Accel_time_n[iVar] = val_fea[iVar+2*nVar];
//    }
//  }
//
//  if (fsi_analysis) {
//    Solution_Pred       =  new su2double [nVar];
//    Solution_Pred_Old   =  new su2double [nVar];
//    Solution_BGS_k      =  new su2double [nVar];
//    for (iVar = 0; iVar < nVar; iVar++) {
//      Solution_Pred[iVar]     = val_fea[iVar];
//      Solution_Pred_Old[iVar] = val_fea[iVar];
//      Solution_BGS_k[iVar]    = 0.0;
//    }
//  }
//
//  /*--- This variable is not "ours", careful not to leak memory ---*/
//  if (Solution_Old == NULL)
//  {
//    /*--- If we are going to use incremental analysis, we need a way to store the old solution ---*/
//    if (incremental_load && nonlinear_analysis) {
//      Solution_Old = new su2double [nVar];
//      for (iVar = 0; iVar < nVar; iVar++) Solution_Old[iVar] = 0.0;
//    }
//    /*--- If we are running a discrete adjoint iteration, we need this vector for cross-dependencies ---*/
//    else if (discrete_adjoint && fsi_analysis) {
//      Solution_Old = new su2double [nVar];
//      for (iVar = 0; iVar < nVar; iVar++) Solution_Old[iVar] = val_fea[iVar];
//    }
//  }
//
//  /*--- Body residual ---*/
//  if (body_forces) {
//    Residual_Ext_Body = new su2double [nVar];
//    for (iVar = 0; iVar < nVar; iVar++) Residual_Ext_Body[iVar] = 0.0;
//  }
//
//  if (refgeom) Reference_Geometry = new su2double [nVar];
//
//  if (prestretch_fem)  Prestretch = new su2double [nVar];
//
//}

CFEAVariable::CFEAVariable(Idx_t npoint, Idx_t ndim, Idx_t nvar, CConfig *config)
  : CVariable(npoint, ndim, nvar, config) {

  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);
  bool body_forces        = config->GetDeadLoad();
  bool incremental_load   = config->GetIncrementalLoad();
  bool prestretch_fem     = config->GetPrestretch();  // Structure is prestretched
  bool discrete_adjoint   = config->GetDiscrete_Adjoint();
  bool refgeom            = config->GetRefGeom(); // Reference geometry needs to be stored
  bool dynamic_analysis   = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fsi_analysis       = config->GetFSI_Simulation();

  VonMises_Stress.resize(nPoint) = 0.0;

  if (nDim==2) Stress.resize(nPoint,3);
  else         Stress.resize(nPoint,6);

  /*--- Initialization of variables ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    Solution[iVar] = val_fea[iVar];
//  }

  if (dynamic_analysis) {
    Solution_Vel.resize(nPoint,nVar);
    Solution_Accel.resize(nPoint,nVar);
    Solution_Vel_time_n.resize(nPoint,nVar);
    Solution_Accel_time_n.resize(nPoint,nVar);
//    for (iVar = 0; iVar < nVar; iVar++) {
//      Solution_Vel[iVar]          = val_fea[iVar+nVar];
//      Solution_Vel_time_n[iVar]   = val_fea[iVar+nVar];
//      Solution_Accel[iVar]        = val_fea[iVar+2*nVar];
//      Solution_Accel_time_n[iVar] = val_fea[iVar+2*nVar];
//    }
  }

  if (fsi_analysis) {
    Solution_Pred.resize(nPoint,nVar);
    Solution_Pred_Old.resize(nPoint,nVar);
    Solution_BGS_k.resize(nPoint,nVar) = 0.0;
//    for (iVar = 0; iVar < nVar; iVar++) {
//      Solution_Pred[iVar]     = val_fea[iVar];
//      Solution_Pred_Old[iVar] = val_fea[iVar];
//    }
  }

  /*--- If we are going to use incremental analysis, we need a way to store the old solution ---*/
  if (incremental_load && nonlinear_analysis) {
    Solution_Old.resize(nPoint,nVar) = 0.0;
  }
  /*--- If we are running a discrete adjoint iteration, we need this vector for cross-dependencies ---*/
  else if (discrete_adjoint && fsi_analysis) {
    Solution_Old.resize(nPoint,nVar);
//    for (iVar = 0; iVar < nVar; iVar++) Solution_Old[iVar] = val_fea[iVar];
  }

  /*--- Body residual ---*/
  if (body_forces) Residual_Ext_Body.resize(nPoint,nVar) = 0.0;

  if (refgeom) Reference_Geometry.resize(nPoint,nVar);

  if (prestretch_fem)  Prestretch.resize(nPoint,nVar);

}
