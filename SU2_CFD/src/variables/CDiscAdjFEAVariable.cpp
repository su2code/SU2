/*!
 * \file CDiscAdjFEAVariable.cpp
 * \brief Definition of the variables for FEM adjoint elastic structural problems.
 * \author Ruben Sanchez
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

#include "../../include/variables/CDiscAdjFEAVariable.hpp"


CDiscAdjFEAVariable::CDiscAdjFEAVariable(const su2double *disp, const su2double *vel, const su2double *accel, unsigned long npoint,
  unsigned long ndim, unsigned long nvar, bool unsteady, CConfig *config) : CVariable(npoint, ndim, nvar, config) {

  bool fsi = config->GetFSI_Simulation();

  Solution_Direct.resize(nPoint,nVar);

  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = disp[iVar];

  if (fsi) {
    Cross_Term_Derivative.resize(nPoint,nDim) = su2double(0.0);
    Geometry_CrossTerm_Derivative.resize(nPoint,nDim) = su2double(0.0);
    
    Solution_BGS.resize(nPoint,nDim) = su2double(0.0);
  }
  
  if (config->GetMultizone_Problem())
    Solution_BGS_k.resize(nPoint,nDim) = su2double(0.0);

  /*--- Nothing else to allocate ---*/
  if (!unsteady) return;


  Dynamic_Derivative.resize(nPoint,nVar) = su2double(0.0);
  Dynamic_Derivative_n.resize(nPoint,nVar) = su2double(0.0);
  Dynamic_Derivative_Vel.resize(nPoint,nVar) = su2double(0.0);
  Dynamic_Derivative_Vel_n.resize(nPoint,nVar) = su2double(0.0);
  Dynamic_Derivative_Accel.resize(nPoint,nVar) = su2double(0.0);
  Dynamic_Derivative_Accel_n.resize(nPoint,nVar) = su2double(0.0);

  Solution_Direct_Vel.resize(nPoint,nVar) = su2double(0.0);
  Solution_Direct_Accel.resize(nPoint,nVar) = su2double(0.0);

  Solution_Vel.resize(nPoint,nVar);
  Solution_Accel.resize(nPoint,nVar);

  Solution_Old_Vel.resize(nPoint,nVar) = su2double(0.0);
  Solution_Old_Accel.resize(nPoint,nVar) = su2double(0.0);

  Solution_Vel_time_n.resize(nPoint,nVar) = su2double(0.0);
  Solution_Accel_time_n.resize(nPoint,nVar) = su2double(0.0);

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) {
      Solution_Vel(iPoint,iVar) = vel[iVar];
      Solution_Accel(iPoint,iVar) = accel[iVar];
    }
  }

}

void CDiscAdjFEAVariable::Set_OldSolution_Vel() { Solution_Old_Vel = Solution_Vel; }

void CDiscAdjFEAVariable::Set_OldSolution_Accel() { Solution_Old_Accel = Solution_Accel; }
