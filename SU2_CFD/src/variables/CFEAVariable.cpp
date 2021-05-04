/*!
 * \file CFEAVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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


CFEAVariable::CFEAVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CVariable(npoint, ndim, nvar, config) {

  bool nonlinear_analysis = (config->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE);
  bool body_forces        = config->GetDeadLoad();
  bool incremental_load   = config->GetIncrementalLoad();
  bool prestretch_fem     = config->GetPrestretch();  // Structure is prestretched
  bool discrete_adjoint   = config->GetDiscrete_Adjoint();
  bool refgeom            = config->GetRefGeom(); // Reference geometry needs to be stored
  bool dynamic_analysis   = config->GetTime_Domain();
  bool multizone          = config->GetMultizone_Problem();
  bool fsi_analysis       = config->GetFSI_Simulation() || multizone;

  VonMises_Stress.resize(nPoint) = su2double(0.0);

  if (nDim==2) Stress.resize(nPoint,3);
  else         Stress.resize(nPoint,6);

  /*--- Initialization of variables ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = val_fea[iVar];

  if (dynamic_analysis) {
    Solution_Vel.resize(nPoint,nVar);
    Solution_Accel.resize(nPoint,nVar);

    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      for (unsigned long iVar = 0; iVar < nVar; iVar++) {
        Solution_Vel(iPoint,iVar) = val_fea[iVar+nVar];
        Solution_Accel(iPoint,iVar) = val_fea[iVar+2*nVar];
      }
    }
    Solution_Vel_time_n = Solution_Vel;
    Solution_Accel_time_n = Solution_Accel;

    if(config->GetMultizone_Problem() && config->GetAD_Mode()) {
      AD_Vel_InputIndex.resize(nPoint,nDim) = -1;
      AD_Vel_OutputIndex.resize(nPoint,nDim) = -1;
      AD_Vel_Time_n_InputIndex.resize(nPoint,nDim) = -1;
      AD_Vel_Time_n_OutputIndex.resize(nPoint,nDim) = -1;
      AD_Accel_InputIndex.resize(nPoint,nDim) = -1;
      AD_Accel_OutputIndex.resize(nPoint,nDim) = -1;
      AD_Accel_Time_n_InputIndex.resize(nPoint,nDim) = -1;
      AD_Accel_Time_n_OutputIndex.resize(nPoint,nDim) = -1;
    }
  }

  if (fsi_analysis) {
    Solution_Pred = Solution;
    Solution_Pred_Old = Solution;
    if (dynamic_analysis) Solution_Vel_Pred = Solution_Vel;
  }

  /*--- If we are going to use incremental analysis, we need a way to store the old solution ---*/

  if (incremental_load && nonlinear_analysis) Solution_Old.resize(nPoint,nVar) = su2double(0.0);

  /*--- If we are running a discrete adjoint iteration, we need this vector for cross-dependencies ---*/

  else if (discrete_adjoint && fsi_analysis) Solution_Old = Solution;

  /*--- Body residual ---*/
  if (body_forces) Residual_Ext_Body.resize(nPoint,nVar) = su2double(0.0);

  if (refgeom) Reference_Geometry.resize(nPoint,nVar);

  if (prestretch_fem) Prestretch.resize(nPoint,nVar);

  if (multizone) Set_BGSSolution_k();

  if (config->GetTopology_Optimization()) {
    nAuxVar = 1;
    AuxVar.resize(nPoint);
  }
}

void CFEAVariable::SetSolution_Vel_time_n() { Solution_Vel_time_n = Solution_Vel; }

void CFEAVariable::SetSolution_Accel_time_n() { Solution_Accel_time_n = Solution_Accel; }

void CFEAVariable::Register_femSolution(bool input, bool push_index, bool dynamic) {
  SU2_OMP_FOR_STAT(roundUpDiv(nPoint,omp_get_num_threads()))
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for(unsigned long iVar=0; iVar<nDim; ++iVar) {
      if(input) {
        if(push_index) {
          AD::RegisterInput(Solution(iPoint,iVar));
          if (dynamic) {
            AD::RegisterInput(Solution_Vel(iPoint,iVar));
            AD::RegisterInput(Solution_Accel(iPoint,iVar));
          }
        }
        else {
          AD::RegisterInput(Solution(iPoint,iVar), false);
          AD::SetIndex(AD_InputIndex(iPoint,iVar), Solution(iPoint,iVar));
          if (dynamic) {
            AD::RegisterInput(Solution_Vel(iPoint,iVar), false);
            AD::SetIndex(AD_Vel_InputIndex(iPoint,iVar), Solution_Vel(iPoint,iVar));
            AD::RegisterInput(Solution_Accel(iPoint,iVar), false);
            AD::SetIndex(AD_Accel_InputIndex(iPoint,iVar), Solution_Accel(iPoint,iVar));
          }
        }
      }
      else {
        AD::RegisterOutput(Solution(iPoint,iVar));
        if(!push_index)
          AD::SetIndex(AD_OutputIndex(iPoint,iVar), Solution(iPoint,iVar));
        if (dynamic) {
          AD::RegisterOutput(Solution_Vel(iPoint,iVar));
          if(!push_index)
            AD::SetIndex(AD_Vel_OutputIndex(iPoint,iVar), Solution_Vel(iPoint,iVar));
          AD::RegisterOutput(Solution_Accel(iPoint,iVar));
          if(!push_index)
            AD::SetIndex(AD_Accel_OutputIndex(iPoint,iVar), Solution_Accel(iPoint,iVar));
        }
      }
    }
  }
  END_SU2_OMP_FOR
}

void CFEAVariable::Register_femSolution_time_n(bool input, bool push_index, bool dynamic) {
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for(unsigned long iVar=0; iVar<nVar; ++iVar) {
      if(input) {
        if(push_index) {
          AD::RegisterInput(Solution_time_n(iPoint,iVar));
          if (dynamic) {
            AD::RegisterInput(Solution_Vel_time_n(iPoint,iVar));
            AD::RegisterInput(Solution_Accel_time_n(iPoint,iVar));
          }
        }
        else {
          AD::RegisterInput(Solution_time_n(iPoint,iVar), false);
          AD::SetIndex(AD_Time_n_InputIndex(iPoint,iVar), Solution_time_n(iPoint,iVar));
          if (dynamic) {
            AD::RegisterInput(Solution_Vel_time_n(iPoint,iVar), false);
            AD::SetIndex(AD_Vel_Time_n_InputIndex(iPoint,iVar), Solution_Vel_time_n(iPoint,iVar));
            AD::RegisterInput(Solution_Accel_time_n(iPoint,iVar), false);
            AD::SetIndex(AD_Accel_Time_n_InputIndex(iPoint,iVar), Solution_Accel_time_n(iPoint,iVar));
          }
        }
      }
      else {
        AD::RegisterOutput(Solution_time_n(iPoint,iVar));
        if(!push_index)
          AD::SetIndex(AD_Time_n_OutputIndex(iPoint,iVar), Solution_time_n(iPoint,iVar));
        if (dynamic) {
          AD::RegisterOutput(Solution_Vel_time_n(iPoint,iVar));
          if(!push_index)
            AD::SetIndex(AD_Vel_Time_n_OutputIndex(iPoint,iVar), Solution_Vel_time_n(iPoint,iVar));
          AD::RegisterOutput(Solution_Accel_time_n(iPoint,iVar));
          if(!push_index)
            AD::SetIndex(AD_Accel_Time_n_OutputIndex(iPoint,iVar), Solution_Accel_time_n(iPoint,iVar));
        }
      }
    }
  }
}

void CFEAVariable::SetfemAdjointSolution(unsigned long iPoint, const su2double *adj_sol, bool dynamic) {
  if (!dynamic) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }
  else {
    for (unsigned long iVar = 0; iVar < 3*nDim; iVar++)
      SU2_TYPE::SetDerivative(Solution(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }
}

void CFEAVariable::SetfemAdjointSolution_LocalIndex(unsigned long iPoint, const su2double *adj_sol, bool dynamic) {
  for (unsigned long iVar = 0; iVar < nDim; iVar++) {
    AD::SetDerivative(AD_OutputIndex(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
    if (dynamic) {
      AD::SetDerivative(AD_Vel_OutputIndex(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar+nDim]));
      AD::SetDerivative(AD_Accel_OutputIndex(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar+2*nDim]));
    }
  }
}

void CFEAVariable::GetfemAdjointSolution(unsigned long iPoint, su2double *adj_sol, bool dynamic) {
  for (unsigned long iVar = 0; iVar < nDim; iVar++) {
    adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution(iPoint,iVar));
    if (dynamic) {
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Vel(iPoint,iVar));
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Accel(iPoint,iVar));
    }
  }
}

void CFEAVariable::GetfemAdjointSolution_LocalIndex(unsigned long iPoint, su2double *adj_sol, bool dynamic) {
  for (unsigned long iVar = 0; iVar < nVar; iVar++) {
    adj_sol[iVar] = AD::GetDerivative(AD_InputIndex(iPoint,iVar));
    if (dynamic) {
      adj_sol[iVar] = AD::GetDerivative(AD_Vel_InputIndex(iPoint,iVar));
      adj_sol[iVar] = AD::GetDerivative(AD_Accel_InputIndex(iPoint,iVar));
    }
  }
}

void CFEAVariable::GetfemAdjointSolution_time_n(unsigned long iPoint, su2double *adj_sol, bool dynamic) {
  for (unsigned long iVar = 0; iVar < nVar; iVar++) {
    adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n(iPoint,iVar));
    if (dynamic) {
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Vel_time_n(iPoint,iVar));
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Accel_time_n(iPoint,iVar));
    }
  }
}

void CFEAVariable::GetfemAdjointSolution_time_n_LocalIndex(unsigned long iPoint, su2double *adj_sol, bool dynamic) {
  for (unsigned long iVar = 0; iVar < nDim; iVar++) {
    adj_sol[iVar] = AD::GetDerivative(AD_Time_n_InputIndex(iPoint,iVar));
    if (dynamic) {
      adj_sol[iVar] = AD::GetDerivative(AD_Vel_Time_n_InputIndex(iPoint,iVar));
      adj_sol[iVar] = AD::GetDerivative(AD_Accel_Time_n_InputIndex(iPoint,iVar));
    }
  }
}
