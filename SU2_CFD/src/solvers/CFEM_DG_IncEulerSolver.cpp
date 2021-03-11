/*!
 * \file CFEM_DG_IncEulerSolver.cpp
 * \brief Main subroutines for solving finite element incompressible Euler flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
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


#include "../../include/solvers/CFEM_DG_IncEulerSolver.hpp"

CFEM_DG_IncEulerSolver::CFEM_DG_IncEulerSolver(CGeometry      *geometry,
                                               CConfig        *config,
                                               unsigned short iMesh)
  : CFEM_DG_SolverBase(geometry, config, iMesh) {

  /*--- Set the number of variables. ---*/
  nVar = nDim + 1;

  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual_RMS.resize(nVar,1.e-35);
  Residual_Max.resize(nVar,1.e-35);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Determine the local and total number of pressure DOFs in the simulation. ---*/
  nPDOFsLocOwned = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    nPDOFsLocOwned += volElem[i].standardElemP->GetNDOFs();

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPDOFsLocOwned, &nPDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#else
  nPDOFsGlobal = nPDOFsLocOwned;
#endif

  /*--- Determine the communication pattern. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  if( !DGGeometry) SU2_MPI::Error(string("Dynamic cast failed"), CURRENT_FUNCTION);
  Prepare_MPI_Communication(DGGeometry, config);


  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

CFEM_DG_IncEulerSolver::~CFEM_DG_IncEulerSolver(void) {
}
