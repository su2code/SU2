/*!
 * \file CVolumeElementFEM_DG.cpp
 * \brief Implementations of the member functions of CVolumeElementFEM_DG.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/geometry/fem_grid/CVolumeElementFEM_DG.hpp"

/*---------------------------------------------------------------------*/
/*---      Public member functions of CVolumeElementFEM_DG.         ---*/
/*---------------------------------------------------------------------*/

void CVolumeElementFEM_DG::AllocateCompressibleFlowVar(CConfig        *config,
                                                       unsigned short nVar){

  /*--- Determine the number of solution DOFs and allocate the
        memory for the regular solution DOFs and the working
        variables of the solution DOFs. ---*/
  const unsigned short nDOFs = standardElemFlow->GetNDOFs();
  solDOFs.resize(nDOFs,nVar);
  solDOFsWork.resize(nDOFs,nVar);

  /*--- Allocate the memory for the transformation matrix between
        conservative and entropy variables in the integration points.
        As this matrix is symmetric, only the upper half (or lower
        half) needs to be stored. ---*/
  const unsigned short nIntPad = standardElemFlow->GetNIntegrationPad();
  dUdVInt.resize(nIntPad,nVar*(nVar+1)/2);

  /*--- Check for the classical Runge-Kutta time integration scheme.
        If used, allocate the memory to store the new state in
        solDOFsAux for an owned element. ---*/
  if(config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT)
    if( elemIsOwned ) solDOFsAux.resize(nDOFs,nVar);

  /*--- Check for the ADER-DG time integration scheme and allocate
        the memory for the additional solution variables. ---*/
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {

    /*--- The predictor solution of all time DOFs. ---*/
    const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();
    solDOFsADERPredictor.resize(nTimeDOFs);
    for(unsigned short i=0; i<nTimeDOFs; ++i)
      solDOFsADERPredictor[i].resize(nDOFs,nVar);

    /*--- Check if the element is adjacent to elements of a lower
          time level and allocate the solution to store this data. ---*/
    if( elemAdjLowerTimeLevel ) solDOFsAux.resize(nDOFs,nVar);
  }
}

void CVolumeElementFEM_DG::AllocateResiduals(CConfig        *config,
                                             unsigned short nVar) {

  /*--- Determine the number of padded solution DOFs and
        allocate the memory for the regular residual. ---*/
  const unsigned short nDOFsPad = standardElemFlow->GetNDOFsPad();
  resDOFs.resize(nDOFsPad,nVar);

  /*--- For ADER another residual must be allocated. ---*/
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG)
    resTotDOFsADER.resize(nDOFsPad,nVar);
}

vector<ColMajorMatrix<su2double> > &CVolumeElementFEM_DG::ComputeGradSolIntPoints(void) {

  /*--- Perform the gemm call to compute the gradient of the solution in
        the integration points and return the vector of matrices. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->GradSolIntPoints(solDOFsWork, standardElemFlow->workGradSolInt[thread]);
  return standardElemFlow->workGradSolInt[thread];
}

ColMajorMatrix<su2double> &CVolumeElementFEM_DG::ComputeSolIntPoints(ColMajorMatrix<su2double> &sol) {

  /*--- Perform the gemm call to compute the solution in the
        integration points and return the matrix. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->SolIntPoints(sol, standardElemFlow->workSolInt[thread]);
  return standardElemFlow->workSolInt[thread];
}

void CVolumeElementFEM_DG::ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt) {
  standardElemFlow->ResidualBasisFunctions(scalarDataInt, resDOFs);
}

void CVolumeElementFEM_DG::ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt){
  standardElemFlow->ResidualGradientBasisFunctions(vectorDataInt, resDOFs);
}

void CVolumeElementFEM_DG::SetConstantSolution(const su2double *sol,
                                               unsigned short  nVar,
                                               unsigned short  startInd) {

  /*--- Easier storage of the number of DOFs. ---*/
  const unsigned short nDOFs = solDOFs.rows();

  /*--- Determine the inverse value of the first basis function
        for this element. ---*/
  const passivedouble invBasis0 = 1.0/standardElemFlow->ValBasis0();

  /*--- Loop over the number of variables to set. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {
    const unsigned short ii = iVar+startInd;

    /*--- Initialize the solution to zero. ---*/
    SU2_OMP_SIMD
    for(unsigned short i=0; i<nDOFs; ++i)
      solDOFs(i,ii) = 0.0;

    /*--- Set the first entry of solDOFs to the constant solution. ---*/
    solDOFs(0,ii) = sol[iVar]*invBasis0;
  }
}
