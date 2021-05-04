/*!
 * \file CVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
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


#include "../../include/variables/CVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"


CVariable::CVariable(unsigned long npoint, unsigned long nvar, CConfig *config) {

  /*--- Initialize the number of solution variables. This version
   of the constructor will be used primarily for converting the
   restart files into solution files (SU2_SOL). ---*/
  nPoint = npoint;
  nVar = nvar;

  /*--- Allocate the solution array. ---*/
  Solution.resize(nPoint,nVar) = su2double(0.0);

  if (config->GetMultizone_Problem())
    Solution_BGS_k.resize(nPoint,nVar) = su2double(0.0);

}

CVariable::CVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config, bool adjoint) {

  /*--- Initializate the number of dimension and number of variables ---*/
  nPoint = npoint;
  nDim = ndim;
  nVar = nvar;

  /*--- Allocate fields common to all problems. Do not allocate fields
   that are specific to one solver, i.e. not common, in this class. ---*/
  Solution.resize(nPoint,nVar) = su2double(0.0);

  Solution_Old.resize(nPoint,nVar) = su2double(0.0);

  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) {
    Solution_time_n.resize(nPoint,nVar);
    Solution_time_n1.resize(nPoint,nVar);
  }
  else if (config->GetTime_Domain()) {
    Solution_time_n.resize(nPoint,nVar) = su2double(0.0);
  }

  if (config->GetMultizone_Problem() && config->GetDiscrete_Adjoint() && !adjoint) {
    AD_InputIndex.resize(nPoint,nVar) = -1;
    AD_OutputIndex.resize(nPoint,nVar) = -1;
    if (config->GetTime_Domain()) {
      AD_Time_n_InputIndex.resize(nPoint,nVar) = -1;
      AD_Time_n1_InputIndex.resize(nPoint,nVar) = -1;
    }
  }

  if (config->GetMultizone_Problem())
    Solution_BGS_k.resize(nPoint,nVar) = su2double(0.0);
}

void CVariable::Set_OldSolution() {
  assert(Solution_Old.size() == Solution.size());
  parallelCopy(Solution.size(), Solution.data(), Solution_Old.data());
}

void CVariable::Set_Solution() {
  assert(Solution.size() == Solution_Old.size());
  parallelCopy(Solution_Old.size(), Solution_Old.data(), Solution.data());
}

void CVariable::Set_Solution_time_n() {
  assert(Solution_time_n.size() == Solution.size());
  parallelCopy(Solution.size(), Solution.data(), Solution_time_n.data());
}

void CVariable::Set_Solution_time_n1() {
  assert(Solution_time_n1.size() == Solution_time_n.size());
  parallelCopy(Solution_time_n.size(), Solution_time_n.data(), Solution_time_n1.data());
}

void CVariable::Set_BGSSolution_k() {
  assert(Solution_BGS_k.size() == Solution.size());
  parallelCopy(Solution.size(), Solution.data(), Solution_BGS_k.data());
}

void CVariable::Restore_BGSSolution_k() {
  assert(Solution.size() == Solution_BGS_k.size());
  parallelCopy(Solution_BGS_k.size(), Solution_BGS_k.data(), Solution.data());
}

void CVariable::SetExternalZero() { parallelSet(External.size(), 0.0, External.data()); }

namespace {
  void RegisterContainer(bool input, bool push_index, su2activematrix& variable, su2matrix<int>& ad_index) {
    const auto nPoint = variable.rows();
    SU2_OMP_FOR_STAT(roundUpDiv(nPoint,omp_get_num_threads()))
    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      for(unsigned long iVar=0; iVar<variable.cols(); ++iVar) {

        if (input) AD::RegisterInput(variable(iPoint,iVar), push_index);
        else AD::RegisterOutput(variable(iPoint,iVar));

        if (!push_index) AD::SetIndex(ad_index(iPoint,iVar), variable(iPoint,iVar));
      }
    }
    END_SU2_OMP_FOR
  }
}

void CVariable::RegisterSolution(bool input, bool push_index) {
  RegisterContainer(input, push_index, Solution, input? AD_InputIndex : AD_OutputIndex);
}

void CVariable::RegisterSolution_time_n(bool push_index) {
  RegisterContainer(true, push_index, Solution_time_n, AD_Time_n_InputIndex);
}

void CVariable::RegisterSolution_time_n1(bool push_index) {
  RegisterContainer(true, push_index, Solution_time_n1, AD_Time_n1_InputIndex);
}
